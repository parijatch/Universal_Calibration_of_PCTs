# tau = 0 ------

plot.df <- run.simu.mvt(n=400, nu=3, rho=0.1, effect.size=0)[,7:10]
colnames(plot.df) <- c("Pareto", "Cauchy", "Cauchy+", "Fréchet")

fig.1 <- splom(~plot.df[,c(4,3,2,1)], 
      # 1. Force the aspect ratio of the panels to be square
      aspect = 1, 
      pch=20,
      col="black",
      cex=0.5,
      # 2. Force the axis limits for every variable to be exactly (0, 1)
      pscales = 0,
      scales = list(draw = FALSE),
      # 3. Your custom panel function with the diagonal line
      panel = function(x, y, ...) {
          panel.splom(x, y, ...)
          panel.abline(a = 0, b = 1, lty = 2)
      }, 
      xlab=TeX("$\\tau = 0$"))

print(fig.1)

# tau = 4 ------

plot.df <- run.simu.mvt(n=400, nu=3, rho=0.1, effect.size=4)[,7:10]
colnames(plot.df) <- c("Pareto", "Cauchy", "Cauchy+", "Fréchet")

fig.2 <- splom(~plot.df[,c(4,3,2,1)], 
      # 1. Force the aspect ratio of the panels to be square
      aspect = 1, 
      pch=20,
      col="black",
      cex=0.5,
      # 2. Force the axis limits for every variable to be exactly (0, 1)
      pscales = 0,
      scales = list(draw = FALSE),
      # 3. Your custom panel function with the diagonal line
      panel = function(x, y, ...) {
          panel.splom(x, y, ...)
          panel.abline(a = 0, b = 1, lty = 2)
      }, 
      xlab=TeX("$\\tau = 4$"))

print(fig.2)