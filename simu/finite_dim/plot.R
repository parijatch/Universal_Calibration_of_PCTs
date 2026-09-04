library(plyr)
library(data.table)
library(lattice)
library(latticeExtra)
library(latex2exp)
library(future)
library(future.apply)
library(progressr)

simu.filenames <- list.files(pattern=".csv.gz", full.names=TRUE)
cat(sprintf("%d simulation files found\n", length(simu.filenames)))
stopifnot(length(simu.filenames) > 0)

# Setup Parameters ------
alpha.power <- 0.05
inv.alpha.vec <- 10**(seq(1, 3, length.out=11))
alpha.vec <- 1 / inv.alpha.vec

# Initialize parallel processing
n_cores <- min(40, max(1, availableCores() - 1))
cat(sprintf("Aggregating data in parallel using %d cores...\n", n_cores))

plan(multisession, workers = n_cores)
handlers(global=TRUE)

with_progress({
  p <- progressor(along = simu.filenames)
  
  summary.list <- future_lapply(simu.filenames, function(.s) {
    .dt <- fread(.s)
    
    .pow.summary <- .dt[, .(
      Pareto = mean(pval.pareto < alpha.power),
      Cauchy = mean(pval.cauchy < alpha.power),
      `Cauchy+` = mean(pval.cauchy.pos < alpha.power),
      Frechet = mean(pval.frechet < alpha.power),
      LRT = mean(pval.LRT < alpha.power),
      n = .N
    ), by = .(mu.type, num.ones, effect.size, d, nu, cor.type, rho)]
    
    .power.dt <- melt(.pow.summary, 
                      measure.vars = c("Pareto", "Cauchy", "Cauchy+", "Frechet", "LRT"),
                      variable.name = "method", value.name = "power")
    
    # -- CALIBRATION SUMMARY--
    .null.dt <- .dt[(mu.type == "A" & effect.size == 0) | (mu.type == "B" & num.ones == 0)]
    .calib.dt <- NULL
    
    if (nrow(.null.dt) > 0) {
      calib.alpha.list <- lapply(alpha.vec, function(.a) {
        .null.summary <- .null.dt[, .(
          Pareto = mean(pval.pareto < .a),
          Cauchy = mean(pval.cauchy < .a),
          `Cauchy+` = mean(pval.cauchy.pos < .a),
          Frechet = mean(pval.frechet < .a),
          LRT = mean(pval.LRT < .a),
          alpha = .a,
          n = .N
        ), by = .(mu.type, num.ones, effect.size, d, nu, cor.type, rho)]
        return(.null.summary)
      })
      
      .calib.wide <- rbindlist(calib.alpha.list)
      .calib.dt <- melt(.calib.wide, 
                        measure.vars = c("Pareto", "Cauchy", "Cauchy+", "Frechet", "LRT"),
                        variable.name = "method", value.name = "rej")
    }
    
    rm(.dt, .null.dt)
    gc()
    
    p() # Update progress bar
    return(list(power = .power.dt, calib = .calib.dt))
  }, future.seed = TRUE)
})

# Combine the aggregated datasets
power.df <- rbindlist(lapply(summary.list, function(x) x$power))
calib.df <- rbindlist(lapply(summary.list, function(x) x$calib))
ratio.power.df <-power.df
rm(summary.list)
gc()

# Format Calibration Data & Plot ------
default_d <- max(as.numeric(as.character(calib.df$d)))

# Isolate just the Type-A null to prevent duplicate points drawing multiple lines
calib.df <- calib.df[method != "LRT" & d == default_d & mu.type == "A"]
calib.df[, cor.type := factor(cor.type)]
calib.df[, d := factor(d)]
calib.df[, nu := factor(nu)]
calib.df[, rho := factor(rho)]
calib.df[, method := factor(method, levels=c("Cauchy+", "Cauchy", "Frechet", "Pareto"))]

log_x_scale <- list(log = 10, equispaced.log = FALSE)
layout_calib <- c(length(unique(calib.df$nu)), length(unique(calib.df$rho)))
colors<-c("#2a9d8f","#1d3557","#f4a261","#e63946")
#shapes<-c(17,16,18, 15)
shapes<-c(2,1,5,0)
theme<- list(
  superpose.line=list(col=colors,lwd=2),
  superpose.symbol = list(col = colors, fill=colors,
                          cex = 0.8, pch = shapes)
)
key.calib <- simpleKey(text = levels(calib.df$method), lines = TRUE, points = TRUE)
key.calib$text$col <- colors     
key.calib$lines$col <- colors    
key.calib$points$col <- colors   
key.calib$points$fill <- colors
key.calib$points$pch <- shapes
key.calib$space <- "top"
key.calib$columns <- length(levels(calib.df$method))

fig.autoreg.calib <- useOuterStrips(
  xyplot(rej/alpha ~ 1/alpha | nu * rho, groups = method,
         data=calib.df[cor.type=="autoreg"],
         layout=layout_calib,
         abline=list(h=1, lty="dashed"), type=c("b","g"),
         par.settings=theme,
         key=key.calib,
         scales=list(x=log_x_scale), 
         xlab=TeX("$1/\\alpha$"), ylab=TeX("Empirical Rejection rate / $\\alpha$")), 
  strip = strip.custom(strip.names = TRUE, var.name = TeX("$\\nu$"), sep = "="),
  strip.left = strip.custom(strip.names = TRUE, var.name = TeX("$\\rho$"), sep = "="))
print(fig.autoreg.calib)

fig.exch.calib <- useOuterStrips(
  xyplot(rej/alpha ~ 1/alpha | nu * rho, groups = method,
         data=calib.df[cor.type=="exch"],
         layout=layout_calib,
         abline=list(h=1, lty="dashed"), type=c("b","g"),
         par.settings=theme,
         key=key.calib,
         scales=list(x=log_x_scale),
         xlab=TeX("$1/\\alpha$"), ylab=TeX("Empirical Rejection rate / $\\alpha$")),
  strip = strip.custom(strip.names = TRUE, var.name = TeX("$\\nu$"), sep = "="),
  strip.left = strip.custom(strip.names = TRUE, var.name = TeX("$\\rho$"), sep = "="))
print(fig.exch.calib)

# Format Power Data & Plot ------
setorder(power.df, method, mu.type, num.ones, effect.size, d, nu, cor.type, rho)

# Calculate Relative Power
lrt.df <- power.df[method == "LRT", .(mu.type, num.ones, effect.size, d, nu, cor.type, rho, lrt_power = power)]
power.df <- merge(power.df[method != "LRT"], lrt.df,
                  by = c("mu.type", "num.ones", "effect.size", "d", "nu", "cor.type", "rho"))

power.df[, rel.power := power / lrt_power]
power.df[, nu := factor(nu)]
power.df[, d := factor(d)]
power.df[, rho := factor(rho)]
power.df[, method := factor(as.character(method), levels=c( "Cauchy+", "Cauchy", "Frechet", "Pareto"))]

default_rho <- sort(unique(power.df$rho))[1]
# 
# # --- Type A Plots  ---
 layout_pow_A <- c(length(unique(power.df$nu)), length(unique(power.df$d)))
 pow_theme<- list(
   superpose.line=list(col=colors,lwd=2),
   superpose.symbol = list(col = colors, fill = colors, 
                           cex = 0.8, pch = shapes)
)
 
 key.pow <- simpleKey(text = levels(power.df$method), lines = TRUE, points = TRUE)
 key.pow$text$col <- colors     
 key.pow$lines$col <- colors    
 key.pow$points$col <- colors   
 key.pow$points$fill <- colors
 key.pow$points$pch <- shapes
 key.pow$space <- "top"
 key.pow$columns <- length(levels(power.df$method))

plot_power_A <- function(cor_val) {
  useOuterStrips(
    xyplot(rel.power ~ effect.size | nu + d,
           groups = method,
           data=power.df[cor.type==cor_val & mu.type=="A" & rho==default_rho & effect.size > 0],
           abline=list(h=1, lty="dashed"),
           layout=layout_pow_A, type=c("b","g"),
           par.settings=pow_theme,
           key=key.pow,
           xlab=TeX("Effect Size ($\\tau$)"), ylab="Power relative to LRT"),
    strip = strip.custom(strip.names = TRUE, var.name = TeX("$\\nu$"), sep = "="),
    strip.left = strip.custom(strip.names = TRUE, var.name = TeX("$d$"), sep = "=")
  )
}

fig.pow.autoreg.A <- plot_power_A("autoreg")
print(fig.pow.autoreg.A)

fig.pow.exch.A <- plot_power_A("exch")
print(fig.pow.exch.A)


# --- Type B Plots ---
max_d <- max(as.numeric(as.character(power.df$d)))

power.df.B <- power.df[mu.type == "B" &
                         rho == default_rho &
                         d == max_d &
                         num.ones > 0 &
                         (num.ones %% 5 == 0)]

layout_pow_B <- c(length(unique(power.df.B$nu)), 1)

plot_power_B <- function(cor_val) {
  xyplot(rel.power ~ num.ones | nu,
         groups = method,
         data=power.df.B[cor.type==cor_val],
         abline=list(h=1, lty="dashed"),
         layout=layout_pow_B, type=c("b","g"),
         par.settings=pow_theme,
         key=key.pow,
         xlab="Number of Non-Zero Entries", ylab="Power relative to LRT",
         strip = strip.custom(strip.names = TRUE, var.name = TeX("$\\nu$"), sep = "="))
}

fig.pow.autoreg.B <- plot_power_B("autoreg")
print(fig.pow.autoreg.B)

fig.pow.exch.B <- plot_power_B("exch")
print(fig.pow.exch.B)

 # =====================================================================
 # Relative Power vs Pareto
 # =====================================================================
 
  power.A.df <- ratio.power.df[mu.type == "A" & effect.size > 0 & method != "LRT"]

 # If you want to drop Frechet, use the line below -- CAUTION!! - Plot parameters and dataframe factors need changing
 # when using the line below
  #power.A.df <- ratio.power.df[mu.type == "A" & effect.size > 0 & method != "LRT" & method != "Frechet"]

 pareto.df <- power.A.df[method == "Pareto", .(effect.size, d, nu, cor.type, rho, pareto_power = power)]
 
 ratio.df <- merge(power.A.df[method != "Pareto"], pareto.df, 
                   by = c("effect.size", "d", "nu", "cor.type", "rho"))
 
 ratio.df[, rel.power.pareto := power / pareto_power]
 
 ratio.df[, method := factor(as.character(method), levels=c("Cauchy+", "Cauchy","Frechet"))]
 ratio.df[, nu :=factor(nu)]
 ratio.df[,d:=factor(d)]
 key.ratio <- simpleKey(text = levels(ratio.df$method), lines = TRUE, points = TRUE)
 key.ratio$text$col <- colors[1:3]     
 key.ratio$lines$col <- colors[1:3]   
 key.ratio$points$col <- colors [1:3] 
 key.ratio$points$fill <- colors[1:3]
 key.ratio$points$pch <- shapes[1:3]
 key.ratio$space <- "top"
 key.ratio$columns <- length(levels(ratio.df$method))
 ratio_theme <- list(
   superpose.line = list(col = colors[1:3], lwd = 2),
   superpose.symbol=list(col=colors[1:3],fill=colors[1:3],cex=0.8, pch = shapes[1:3])
 )
 
 # Plotting Function for the Ratio
 plot_ratio_A <- function(cor_val) {
   useOuterStrips(
     xyplot(rel.power.pareto ~ effect.size | nu + d, 
            groups = method, 
            data=ratio.df[cor.type==cor_val & rho==default_rho], 
            abline=list(h=1, lty="dashed"),
            layout=layout_pow_A, type=c("b","g"), 
            par.settings = ratio_theme,
            key=key.ratio,
            xlab=TeX("Effect Size ($\\tau$)"), ylab="Power relative to Pareto",
            ),
     strip = strip.custom(strip.names = TRUE, var.name = TeX("$\\nu$"), sep = "="),
     strip.left = strip.custom(strip.names = TRUE, var.name = TeX("$d$"), sep = "=")
   )
 }
 
 fig.ratio.autoreg.A <- plot_ratio_A("autoreg")
 print(fig.ratio.autoreg.A)
 
 fig.ratio.exch.A <- plot_ratio_A("exch")
 print(fig.ratio.exch.A)
 