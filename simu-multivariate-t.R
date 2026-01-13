library(plyr)
library(MASS)
library(data.table)
library(latex2exp)

source("utils.R")

get.cor <- function(d, rho=0.5, cor.type="exch") {
  stopifnot(cor.type %in% c("exch", "autoreg", "diag"))
  if (cor.type=="exch") {
    cor.mat <- matrix(rep(rho, d*d), nrow=d)
    diag(cor.mat) <- 1
  }
  else if (cor.type=="autoreg") {
    cor.mat <- rho^(abs(outer(1:d, 1:d, FUN="-")))
  } else if (cor.type=="diag") {
    cor.mat <- diag(d)
  }
  return(cor.mat)
} 
  

simu.mvt <- function(d=4, n=1e5, nu=1, cor.mat=diag(d), 
                     effect.size=0, d.nonzero=1) {
  stopifnot(all(diag(cor.mat)==1))
  # multivariate t with cov=cor.mat
  X <- mvrnorm(n=n, mu=rep(0, d), Sigma=cor.mat)
  denom <- replicate(n, sqrt(mean(rnorm(nu)^2)))
  X <- X / denom
  # add signal to the mean
  stopifnot(d.nonzero <= d & d.nonzero > 0)
  mu <- rep(0, d)
  mu[1:d.nonzero] <- 1
  mu <- mu / sqrt(sum(mu^2)) * effect.size
  mu <- t(replicate(n, mu))
  X <- X + mu
  # return pval
  pval <- 1 - pt(X, df=nu)
  return(pval)
}

test.LRT <- function(pval, pval.null, 
                     effect.size, d.nonzero, nu, cor.mat) {
  X <- qt(1 - pval, df=nu)
  X0 <- qt(1 - pval.null, df=nu)
  d <- dim(X)[2]
  n <- dim(X)[1]
  mu <- rep(0, d)
  mu[1:d.nonzero] <- 1
  if (effect.size == 0) {
    # exclude effect.size=0 here because LRT = constant 1 under the null 
    effect.size <- 1e-7
  }
  mu <- mu / sqrt(sum(mu^2)) * effect.size
  mu <- t(replicate(n, mu))
  inv.cor.mat <- solve(cor.mat)
  
  l.fun <- function(.X) {
    apply(.X, 1, function(.x) 1 + .x %*% inv.cor.mat %*% .x / nu)
  }
  
  # for X
  LR.X <- l.fun(X) / l.fun(X - mu)
  # for X0
  LR.X0 <- l.fun(X0) / l.fun(X0-mu)
  # get pvals
  F0 <- ecdf(LR.X0)
  pvals <- 1 - F0(LR.X)
  return(pvals)
}

# Simu ------
run.simu.mvt <- function(n=1e4, d=10, nu=1, d.nonzero=5, cor.type="autoreg", 
                         effect.size.vec=seq(0, 40, length.out=20),
                         alpha=0.05, .plot=TRUE) {
  power.df <- ldply(effect.size.vec, function(effect.size) {
    cor.mat <- get.cor(d, cor.type=cor.type)
    # pvals under effect.size
    X <- simu.mvt(d=d,
                  n=n,
                  cor.mat = cor.mat, 
                  d.nonzero = d.nonzero,
                  effect.size = effect.size, nu=nu)
    # pvals under the null
    X0 <- simu.mvt(d=d,
                   n=n,
                   cor.mat = cor.mat, 
                   d.nonzero = d.nonzero,
                   effect.size = 0, nu=nu)
    pval.pareto <- combine.test(X, method="Pareto")
    pval.cauchy <- combine.test(X, method="Cauchy")
    pval.cauchy.pos <- combine.test(X, method="Cauchy+")
    pval.LRT <- test.LRT(X, X0, effect.size, d.nonzero, nu, cor.mat)
    .df <- data.frame(tau=rep(effect.size, 4),
                      power=c(mean(pval.pareto < alpha), 
                              mean(pval.cauchy < alpha), 
                              mean(pval.cauchy.pos < alpha),
                              mean(pval.LRT < alpha)), 
                      method=c("Pareto", "Cauchy", "Cauchy+","LRT"))
  }, .progress = "text")
  power.df <- data.table(power.df)
  
  if (.plot) {
    power.df[method=="Pareto", plot(tau, power, type="b", ylim=c(0,1), 
                                   xlab=TeX(r"(effect size = $\|\mu\|_2$)"), col="red", 
                                   main=TeX(sprintf(r"(Multivariate-$t_\nu$($\mu$, $\Sigma$) with %s $\Sigma$, $\nu$=%d: dim($\mu$)=%d, $\|\mu\|_0$=%d)", 
                                                    cor.type, nu, d, d.nonzero)), 
                                   sub=TeX(sprintf(r"($\alpha$=%g)", alpha)))]
    power.df[method=="Cauchy", points(tau, power, type="b", col="blue")]
    power.df[method=="Cauchy+", points(tau, power, type="b", col="orange")]
    power.df[method=="LRT", points(tau, power, type="b", col="black")]
    abline(h=alpha, lty=2)
    legend("bottomright", legend=c("Pareto", "Cauchy", "Cauchy+","LRT"), pch=c(1,1), col=c("red", "blue","orange","black"))
  }
  return(power.df)
}

# examples ----
#run.simu.mvt(, .plot=TRUE)

#run.simu.mvt(nu=2, .plot=TRUE)

#run.simu.mvt(d=50, d.nonzero=40, nu=3, .plot=TRUE)

#Relative Power
nu.vec <- c(1, 5, 20)                
d.nonzero.vec <- c(1, 4, 8)         

nrows <- length(nu.vec)
ncols <- length(d.nonzero.vec)
cor.type="autoreg"
par(mfrow = c(nrows, ncols + 1), mar = c(4, 4, 3, 1))

plot.idx <- 1
for (nu in nu.vec) {
  for (d.nonzero in d.nonzero.vec) {
    message(sprintf("Running nu=%d, d.nonzero=%d", nu, d.nonzero))
    power.df <- run.simu.mvt(n=1e4, d=10, nu=nu, d.nonzero=d.nonzero,
                             cor.type=cor.type,
                             effect.size.vec=seq(0, 40, length.out=20),
                             alpha=0.05, .plot=FALSE)
    
    # merge power by method and effect size
    power.wide <- dcast(power.df, tau ~ method, value.var="power")
    
    rel.pareto <- power.wide$Pareto / power.wide$LRT
    rel.cauchy <- power.wide$Cauchy / power.wide$LRT
    plot(power.wide$tau, rel.pareto, type="b", ylim=c(0, max(1.5, rel.pareto, rel.cauchy, na.rm=TRUE)),
         xlab=TeX(r"(effect size = $\|\mu\|_2$)"), ylab="Relative power vs LRT",
         col="red",
         main=TeX(sprintf(r"($\nu$=%d, $\|\mu\|_0$=%d)", nu, d.nonzero)))
    points(power.wide$tau, rel.cauchy, type="b", col="blue")
    abline(h=1, lty=2)
    
    plot.idx <- plot.idx + 1
  }
  plot.new()  # empty plot panel for legend
  legend("center", legend=c("Pareto/LRT", "Cauchy/LRT"), col=c("red", "blue"), pch=1, cex=1.2, bty="n")
}
mtext(
  TeX(sprintf("Relative Power comparison of PCT & CCT for Multivariate-$t_\\nu(\\mu, \\Sigma)$ with %s $\\Sigma$", cor.type)),
  outer = TRUE, cex = 1.5, line = 1
)
par(mfrow=c(1,1))  
