library(plyr)
library(data.table)
library(lattice)
library(latticeExtra)
library(latex2exp)

simu.filenames <- list.files(pattern=".csv.gz", full.names=TRUE)

cat(sprintf("%d simulation files found\n", length(simu.filenames)))

stopifnot(length(simu.filenames) > 0)

# calibration ------
alpha.vec <- 1/c(20, 50, 100, 200, 300, 400, 500)

pvals.df <- rbindlist(llply(simu.filenames, function(.s) {
    .dt <- fread(.s)
    if (!.dt[1, effect.size==0 & d==10]) {
        return(NULL)
    } else {
        .dt <- melt(.dt, 
             measure.vars = c("pval.pareto", 
                              "pval.cauchy", 
                              "pval.cauchy.pos",
                              "pval.frechet", 
                              "pval.LRT"), 
             variable.name = "method", value.name = "pval")
        .dt[, method:=revalue(method, c("pval.pareto"="Pareto", 
                                     "pval.cauchy"="Cauchy", 
                                     "pval.cauchy.pos"="Cauchy+", 
                                     "pval.frechet"="Frechet", 
                                     "pval.LRT"="LRT"))]
        return(.dt)
    }
}, .progress = "text"))

calib.df <- pvals.df[effect.size==0 & method!="LRT" & d==10, 
                     .(rej=rowMeans(sapply(pval, function(.x) .x < alpha.vec)), 
                         alpha=alpha.vec, 
                       n=.N), 
                     by=.(cor.type, nu, rho, method, d)]
rm(pvals.df)

calib.df[, cor.type:=factor(cor.type)]
calib.df[, d:=factor(d)]
calib.df[, nu:=factor(nu)]
calib.df[, rho:=factor(rho)]
calib.df[, method:=factor(method, 
                          levels=c("Pareto", "Cauchy", "Cauchy+", "Frechet"))]


print(summary(calib.df))

fig.autoreg <- useOuterStrips(
    xyplot(rej/alpha ~ 1/alpha | nu * rho, groups = method,
           data=calib.df[cor.type=="autoreg"],
           auto.key = TRUE,
           layout=c(4,2),
           abline=list(h=1, lty="dashed"),
           type=c("b","g"),
           xlab=TeX("$1/\\alpha$"),
           ylab=TeX("Rejection rate / $\\alpha$")), 
    strip = strip.custom(strip.names = TRUE, 
                         var.name = TeX("$\\nu$"), 
                         sep = "="),
    strip.left = strip.custom(strip.names = TRUE, 
                              var.name = TeX("$\\rho$"),
                              sep = "="))
print(fig.autoreg)

fig.exch <- useOuterStrips(
    xyplot(rej/alpha ~ 1/alpha | nu * rho, groups = method,
           data=calib.df[cor.type=="exch"],
           auto.key = TRUE,
           layout=c(4,2),
           abline=list(h=1, lty="dashed"),
           type=c("b","g"),
           xlab=TeX("$1/\\alpha$"),
           ylab=TeX("Rejection rate / $\\alpha$")), 
    strip = strip.custom(strip.names = TRUE, 
                         var.name = TeX("$\\nu$"), 
                         sep = "="),
    strip.left = strip.custom(strip.names = TRUE, 
                              var.name = TeX("$\\rho$"),
                              sep = "="))
print(fig.exch)

# power -------

alpha <- 0.05

power.df <- rbindlist(llply(simu.filenames, function(.s) {
    .dt <- fread(.s)
    .dt <- melt(.dt, 
                measure.vars = c("pval.pareto", 
                                 "pval.cauchy", 
                                 "pval.cauchy.pos",
                                 "pval.frechet", 
                                 "pval.LRT"), 
                variable.name = "method", value.name = "pval")
    .dt[, method:=revalue(method, c("pval.pareto"="Pareto", 
                                    "pval.cauchy"="Cauchy", 
                                    "pval.cauchy.pos"="Cauchy+", 
                                    "pval.frechet"="Frechet", 
                                    "pval.LRT"="LRT"))]
    .power.dt <- .dt[, .(power=mean(pval < alpha), n=.N), 
                by=.(method, effect.size, d, nu, cor.type, rho)]
    return(.power.dt)
}, .progress = "text"))

setorder(power.df, method, effect.size, d, nu, cor.type, rho)

power.df[, nu:=factor(nu)]
power.df[, d:=factor(d)]
power.df[, rho:=factor(rho)]


# relative power

.LRT.power <- power.df[method=="LRT", power]
power.df[, rel.power:=power / .LRT.power]
stopifnot(all(power.df[method=="LRT", rel.power==1]))
power.df <- power.df[method!="LRT"]
power.df[, method:=factor(as.character(method), 
                             levels=c("Pareto", "Cauchy", "Cauchy+", "Frechet"))]


fig.pow.autoreg <- useOuterStrips(
    xyplot(rel.power ~ effect.size | nu + d, groups = method, 
           data=power.df[cor.type=="autoreg" & method!="LRT" & rho==0.1 & effect.size>0], 
           abline=list(h=1, lty="dashed"),
           auto.key = TRUE,
           layout=c(3,2),
           type=c("b","g"),
           xlab=TeX("$\\tau$"), 
           ylab=TeX("Power relative to LRT")),
    strip = strip.custom(strip.names = TRUE, 
                         var.name = TeX("$\\nu$"), 
                         sep = "="),
    strip.left = strip.custom(strip.names = TRUE, 
                              var.name = TeX("$d$"), 
                              sep = "="))

print(fig.pow.autoreg)


fig.pow.exch <- useOuterStrips(
    xyplot(rel.power ~ effect.size | nu + d, groups = method, 
           data=power.df[cor.type=="exch" & method!="LRT" & rho==0.1 & effect.size>0], 
           abline=list(h=1, lty="dashed"),
           auto.key = TRUE,
           layout=c(3,2),
           type=c("b","g"),
           xlab=TeX("$\\tau$"), 
           ylab=TeX("Power relative to LRT")),
    strip = strip.custom(strip.names = TRUE, 
                         var.name = TeX("$\\nu$"), 
                         sep = "="),
    strip.left = strip.custom(strip.names = TRUE, 
                              var.name = TeX("$d$"), 
                              sep = "="),
)

print(fig.pow.exch)
