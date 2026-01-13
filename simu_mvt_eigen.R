library(plyr)
library(MASS)
library(data.table)
library(latex2exp)

source("~/Downloads/MRV_Cauchy_n_Pareto/Code/utils.R", chdir = TRUE)
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
simu.mvt.eig <- function(d, n=1e5, nu=1, cor.mat, effect.size=0, eig.type=c("top", "bottom")) {
  eig.type <- match.arg(eig.type)
  
  # simulate multivariate t
  X <- mvrnorm(n=n, mu=rep(0, d), Sigma=cor.mat)
  denom <- replicate(n, sqrt(mean(rnorm(nu)^2)))
  X <- X / denom
  
  # direction of mu from eig(Σ⁻¹)
  inv.cor <- solve(cor.mat)
  eig <- eigen(inv.cor)
  
  mu.dir <- if (eig.type == "top") eig$vectors[, 1] else eig$vectors[, d]
  
  mu <- mu.dir / sqrt(sum(mu.dir^2)) * effect.size
  
  mu.rep <- t(replicate(n, mu))
  
  X <- X + mu.rep
  pval <-pmin(1 - pt(X, df=nu),1)
  
  
  return(pval)
}

test.LRT.mod <- function(pval, pval.null, effect.size, nu, cor.mat, eig.type) {
  X <- qt(1 - pval, df=nu)
  X0 <- qt(1 - pval.null, df=nu)
  d <- dim(X)[2]
  n <- dim(X)[1]
  inv.cor.mat <- solve(cor.mat)
  
  eig <- eigen(inv.cor.mat)
  mu.dir <- if (eig.type == "top") eig$vectors[, 1] else eig$vectors[, d]
  if (effect.size == 0) effect.size <- 1e-7
  mu <- mu.dir / sqrt(sum(mu.dir^2)) * effect.size
  mu.rep <- matrix(rep(mu, each=n), nrow=n)
  
  quadform <- function(mat) rowSums((mat %*% inv.cor.mat) * mat)
  
  num  <- 1 + quadform(X) / nu
  den  <- 1 + quadform(X - mu.rep) / nu
  num0 <- 1 + quadform(X0) / nu
  den0 <- 1 + quadform(X0 - mu.rep) / nu
  
  LR.X  <- num / den
  LR.X0 <- num0 / den0
  
  # Replace non-finite values with a large finite fallback
  safe_max <- function(x) {
    if (all(!is.finite(x))) return(1e6)
    return(max(x[is.finite(x)], na.rm = TRUE))
  }
  
  LR.X[!is.finite(LR.X)]  <- safe_max(LR.X)
  LR.X0[!is.finite(LR.X0)] <- safe_max(LR.X0)
  
  F0 <- ecdf(LR.X0)
  pvals <- 1 - F0(LR.X)
  
  # Final cleanup
  pvals[is.na(pvals)] <- 1
  return(pvals)
}





# Simulation runner for new μ structure
run.simu.eig <- function(n=1e4, d, nu, effect.size.vec=seq(0, 40, length.out=20),
                         cor.type="autoreg", eig.type=c("top", "bottom"),
                         alpha=0.005, .plot=FALSE) {
  eig.type <- match.arg(eig.type)
  
  power.df <- ldply(effect.size.vec, function(effect.size) {
    cor.mat <- get.cor(d, cor.type=cor.type)
    
    # signal data
    X <- simu.mvt.eig(d=d, n=n, nu=nu, cor.mat=cor.mat,
                      effect.size=effect.size, eig.type=eig.type)
    
    # null data
    X0 <- simu.mvt.eig(d=d, n=n, nu=nu, cor.mat=cor.mat,
                       effect.size=0, eig.type=eig.type)
    
    # compute test pvals
    pval.pareto <- combine.test(X, method="Pareto")
    pval.cauchy <- combine.test(X, method="Cauchy")
    pval.LRT <- test.LRT.mod(X, X0, effect.size=effect.size, nu=nu,
                             cor.mat=cor.mat,eig.type=eig.type)
    if (any(is.na(pval.LRT))) {
      warning(sprintf("NA in LRT p-values for effect.size = %f, nu = %d, d = %d", effect.size, nu, d))
    }
    data.frame(tau=rep(effect.size, 3),
               power=c(mean(pval.pareto < alpha, na.rm=TRUE),
                       mean(pval.cauchy < alpha, na.rm=TRUE),
                       mean(pval.LRT < alpha, na.rm=TRUE)),
               method=c("Pareto", "Cauchy","LRT"))
  })
  
  power.df <- data.table(power.df)
  
  if (.plot) {
    power.df[method=="Pareto", plot(tau, power, type="b", ylim=c(0,1), col="red",
                                    main=TeX(sprintf("μ in %s eig.dir, ν=%d, d=%d", eig.type, nu, d)),
                                    xlab=TeX(r"(effect size = $\|\mu\|_2$)"), ylab="Power")]
    power.df[method=="Cauchy", points(tau, power, type="b", col="blue")]
    power.df[method=="LRT", points(tau, power, type="b", col="black")]
    abline(h=alpha, lty=2)
  }
  
  return(power.df)
}
 ###Plotting###
nu.vec <- c(1, 10, 25)     
d.vec <- c(3, 10, 20)     
eig.type <- "bottom"
cor.type<- "autoreg" 

par(mfrow = c(length(nu.vec), length(d.vec) + 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))

for (nu in nu.vec) {
  for (d in d.vec) {
    message(sprintf("Running nu=%d, d=%d", nu, d))
    
    power.df <- run.simu.eig(n=1e4, d=d, nu=nu, eig.type=eig.type,
                             effect.size.vec=seq(0, 40, length.out=20),
                             cor.type=cor.type, alpha=0.05, .plot=FALSE)
    
    power.wide <- dcast(power.df, tau ~ method, value.var="power")
    
    eps <- 1e-7  # small constant to avoid division by zero
    LRT.safe <- pmax(power.wide$LRT, eps)
    rel.pareto <- power.wide$Pareto / LRT.safe
    rel.cauchy <- power.wide$Cauchy / LRT.safe

    plot(power.wide$tau, rel.pareto, type="b",
         ylim=c(0, max(1.5, rel.pareto, rel.cauchy, na.rm=TRUE)),
         xlab=TeX(r"(effect size = $\|\mu\|_2$)"), ylab="Relative power vs LRT",
         col="red",
         main=TeX(sprintf("ν=%d, d=%d", nu, d)))
    points(power.wide$tau, rel.cauchy, type="b", col="blue")
    abline(h=1, lty=2)
  }
  plot.new()
  legend("center", legend=c("Pareto/LRT", "Cauchy/LRT"), col=c("red", "blue"), pch=1, cex=1.2, bty="n")
}

par(mfrow = c(1, 1))