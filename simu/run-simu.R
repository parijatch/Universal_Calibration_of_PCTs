library(MASS)
library(data.table)
library(future)
library(future.apply)
library(progressr)
plan(multisession)
handlers(global=TRUE)

source("utils.R")

get.cor <- function(d, rho=0.3, cor.type="exch") {
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
                     effect.size=0, eig.type=d) {
  stopifnot(all(diag(cor.mat)==1))
  # multivariate t with cov=cor.mat
  X <- mvrnorm(n=n, mu=rep(0, d), Sigma=cor.mat)
  denom <- sqrt(rchisq(n, df=nu) / nu)
  X <- X / denom
  # add signal to the mean
  stopifnot(1 <= eig.type & eig.type <= d)
  eig.vec <- eigen(cor.mat, symmetric = TRUE)$vectors[,eig.type]
  mu <- effect.size * eig.vec
  X <- X + matrix(rep(mu, n), nrow=n, byrow=TRUE)
  # return two-sided pvals
  pvals <- 2 * pt(abs(X), df=nu, lower.tail = FALSE)
  return(list(pvals=pvals, X=X, effect.size=effect.size,
              mu=mu, nu=nu, cor.mat=cor.mat))
}

test.LRT <- function(X, mu, nu, cor.mat) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  mu.mat <- t(replicate(n, mu))
  X0 <- X - mu.mat
  inv.cor.mat <- solve(cor.mat)
  l.fun <- function(.X) {
    apply(.X, 1, function(.x) 1 + .x %*% inv.cor.mat %*% .x / nu)
  }
  # for X
  LR.X <- l.fun(X) / l.fun(X - mu.mat)
  # for X0
  LR.X0 <- l.fun(X0) / l.fun(X0-mu.mat)
  # get pvals
  F0 <- ecdf(LR.X0)
  pvals <- 1 - F0(LR.X)
  return(pvals)
}

# Simu ------
run.simu.mvt <- function(n=1e4, d=10, nu=1, 
                         cor.type="autoreg", rho=0.5,
                         effect.size=1, eig.type=d) {
  .simu <- simu.mvt(d=d, n=n, nu=nu, 
                    cor.mat=get.cor(d, rho=rho, cor.type=cor.type), 
                    effect.size = effect.size, eig.type=eig.type)
  pval.pareto <- combine.test(.simu$pvals, method="Pareto")
  pval.cauchy <- combine.test(.simu$pvals, method="Cauchy")
  pval.cauchy.pos <- combine.test(.simu$pvals, method="Cauchy+")
  pval.frechet <- combine.test(.simu$pvals, method="Frechet")
  if (effect.size==0) {
    pval.LRT <- runif(n)
  } else {
    pval.LRT <- with(.simu, test.LRT(X, mu, nu, cor.mat))
  }
  data.table(d=d, nu=nu, effect.size=effect.size, 
             cor.type=cor.type, rho=rho, eig.type=eig.type,
             pval.pareto=pval.pareto, 
             pval.cauchy=pval.cauchy,
             pval.cauchy.pos=pval.cauchy.pos, 
             pval.frechet=pval.frechet, 
             pval.LRT=pval.LRT)
}


# config.df <- expand.grid(nu=c(0.5, 1,3, 15), 
#                          d=c(3,10,20),
#                          cor.type=c("exch", "autoreg"), 
#                          rho=c(0.3,0.9),  
#                          effect.size=seq(0, 12, length.out=5))

config.df <- expand.grid(nu=c(3,10,50,1000), 
                         d=c(3,10,20),
                         cor.type=c("exch", "autoreg"), 
                         rho=c(0.1,0.9),  
                         effect.size=seq(0,10,length.out=11))

n <- 1e6
# n <- 1e6

with_progress({
  p <- progressor(steps=nrow(config.df))
  results <- future_lapply(seq_len(nrow(config.df)), function(i) {
    p()
    .config <- config.df[i, ]
    .dt <- with(.config, run.simu.mvt(n=n, d=d, nu=nu, cor.type=cor.type, rho=rho,
                                      effect.size=effect.size, eig.type=d))
    fwrite(.dt, file = tempfile(pattern = "simu-", fileext=".csv.gz", tmpdir = "."))
  }, future.seed=TRUE)
})

cat(sprintf("Saved %d files\n", length(results)))
