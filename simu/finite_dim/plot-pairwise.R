library(lattice)
library(latticeExtra)
library(latex2exp)
library(data.table)

# null ------
plot.df <- run.simu.mvt(n=400, d=100, nu=1, rho=0.2, cor.type="exch", effect.size=0)
plot.df <- plot.df[,c("pval.pareto", "pval.cauchy" ,"pval.cauchy.pos", "pval.frechet")]
colnames(plot.df) <- c("Pareto", "Cauchy", "Cauchy+", "Fréchet")

fig.1 <- splom(~plot.df[,c(4,3,2,1)], 
      # 1. Force the aspect ratio of the panels to be square
      aspect = 1, 
      pch=20,
      col="black",
      cex=0.5,
      # 2. Force the axis limits for every variable to be exactly (0, 1).
      #    prepanel.limits is what sets the limits (it is passed through to
      #    panel.pairs and applied to each variable); pscales only controls
      #    the tick annotation drawn in the diagonal cells.
      prepanel.limits = function(x) c(0, 1),
      pscales = 0,
      scales = list(draw = FALSE),
      # 3. Your custom panel function with the diagonal line
      panel = function(x, y, ...) {
          panel.splom(x, y, ...)
          panel.abline(a = 0, b = 1, lty = 2)
      }, 
      xlab=TeX("$\\tau = 0$"))

print(fig.1)

# alternative ------
plot.df <- run.simu.mvt(n=400, d=100, nu=1, rho=0.2, cor.type="exch", effect.size=20)
plot.df <- plot.df[,c("pval.pareto", "pval.cauchy" ,"pval.cauchy.pos", "pval.frechet")]
colnames(plot.df) <- c("Pareto", "Cauchy", "Cauchy+", "Fréchet")

fig.2 <- splom(~plot.df[,c(4,3,2,1)], 
               # 1. Force the aspect ratio of the panels to be square
               aspect = 1, 
               pch=20,
               col="black",
               cex=0.5,
               # 2. Force the axis limits for every variable to be exactly (0, 1).
               #    prepanel.limits is what sets the limits (it is passed through to
               #    panel.pairs and applied to each variable); pscales only controls
               #    the tick annotation drawn in the diagonal cells.
               prepanel.limits = function(x) c(0, 1),
               pscales = 0,
               scales = list(draw = FALSE),
               # 3. Your custom panel function with the diagonal line
               panel = function(x, y, ...) {
                 panel.splom(x, y, ...)
                 panel.abline(a = 0, b = 1, lty = 2)
               }, 
               xlab=TeX("$\\tau = 20$"))

print(fig.2)
