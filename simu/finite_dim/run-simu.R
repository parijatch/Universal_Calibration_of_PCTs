library(MASS)
library(data.table)
library(future)
library(future.apply)
library(progressr)
library(parallel)

n_cores <- min(40, max(1, availableCores() - 1))
cat(sprintf("Initializing parallel processing with %d cores...\n", n_cores))

plan(multisession, workers = n_cores)
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
                     effect.size=0, mu.type="A", num.ones=0) {
  stopifnot(all(diag(cor.mat)==1))
  # multivariate t with cov=cor.mat
  X <- mvrnorm(n=n, mu=rep(0, d), Sigma=cor.mat)
  denom <- sqrt(rchisq(n, df=nu) / nu)
  X <- X / denom
  
  # Set the location parameter (mu) based on type
  if (mu.type == "A") {
    mu <- (effect.size / sqrt(d)) * rep(1, d)
  } else if (mu.type == "B") {
    stopifnot(num.ones >= 0 && num.ones <= d)
    mu <- c(rep(1, num.ones), rep(0, d - num.ones))
  } else {
    stop("Invalid mu.type. Must be 'A' or 'B'.")
  }
  
  X <- X + matrix(rep(mu, n), nrow=n, byrow=TRUE)
  
  # Return one-sided upper tail pvals
  pvals <- pt(X, df=nu, lower.tail = FALSE)
  return(list(pvals=pvals, X=X, effect.size=effect.size, 
              mu.type=mu.type, num.ones=num.ones, 
              mu=mu, nu=nu, cor.mat=cor.mat))
}

test.LRT <- function(X, mu, nu, cor.mat) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  mu.mat <- matrix(mu, nrow=n, ncol=d, byrow=TRUE)
  
  X0 <- X - mu.mat
  inv.cor.mat <- solve(cor.mat)
  
  l.fun <- function(.X) {
    1 + rowSums((.X %*% inv.cor.mat) * .X) / nu
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
                         effect.size=1, mu.type="A", num.ones=0) {
  .simu <- simu.mvt(d=d, n=n, nu=nu, 
                    cor.mat=get.cor(d, rho=rho, cor.type=cor.type), 
                    effect.size=effect.size, mu.type=mu.type, num.ones=num.ones)
  
  pval.pareto <- combine.test(.simu$pvals, method="Pareto")
  pval.cauchy <- combine.test(.simu$pvals, method="Cauchy")
  pval.cauchy.pos <- combine.test(.simu$pvals, method="Cauchy+")
  pval.frechet <- combine.test(.simu$pvals, method="Frechet")
  
  # Determine if under the global null to calculate LRT
  is.null.hyp <- (mu.type == "A" && effect.size == 0) || (mu.type == "B" && num.ones == 0)
  
  if (is.null.hyp) {
    pval.LRT <- runif(n)
  } else {
    pval.LRT <- with(.simu, test.LRT(X, mu, nu, cor.mat))
  }
  
  data.table(d=d, nu=nu, effect.size=effect.size, 
             mu.type=mu.type, num.ones=num.ones,
             cor.type=cor.type, rho=rho,
             pval.pareto=pval.pareto, 
             pval.cauchy=pval.cauchy,
             pval.cauchy.pos=pval.cauchy.pos, 
             pval.frechet=pval.frechet, 
             pval.LRT=pval.LRT)
}

# Values for the parameter grids
nu_vals <- c(0.5, 1, 10)
d_vals <- c(3, 25, 100)
cor.type_vals <- c("exch", "autoreg")
rho_vals <- c(0.2, 0.8)
effect.size_vals <- seq(0, 50, length.out=11)

# Generate configuration for Type-A parameters
config.df.A <- expand.grid(nu=nu_vals, d=d_vals, cor.type=cor.type_vals, 
                           rho=rho_vals, effect.size=effect.size_vals, 
                           stringsAsFactors=FALSE)
config.df.A$mu.type <- "A"
config.df.A$num.ones <- 0

# Generate configuration for Type-B parameters
grid_base.B <- expand.grid(nu=nu_vals, d=d_vals, cor.type=cor.type_vals, 
                           rho=rho_vals, stringsAsFactors=FALSE)

# For Type-B, map over the base and iterate num.ones from 0 to d for each row
config.df.B <- do.call(rbind, lapply(seq_len(nrow(grid_base.B)), function(i) {
  row <- grid_base.B[i, ]
  do.call(rbind, lapply(0:row$d, function(k) {
    new_row <- row
    new_row$mu.type <- "B"
    new_row$effect.size <- 0  # effect size defaults to 0 for type-B
    new_row$num.ones <- k
    return(new_row)
  }))
}))

# Combine configurations
config.df <- rbind(config.df.A, config.df.B)

n <- 5e5

with_progress({
  p <- progressor(steps=nrow(config.df))
  results <- future_lapply(seq_len(nrow(config.df)), function(i) {
    .config <- config.df[i, ]
    .dt <- with(.config, run.simu.mvt(n=n, d=d, nu=nu, cor.type=cor.type, rho=rho,
                                      effect.size=effect.size, mu.type=mu.type, 
                                      num.ones=num.ones))
    
    fname <- tempfile(pattern = sprintf("simu_type%s_", .config$mu.type), fileext=".csv.gz", tmpdir = ".")
    
    fwrite(.dt, file = fname)
    
    rm(.dt)
    gc()
    p()
    return(fname) 
    
  }, future.seed=TRUE)
})

cat(sprintf("Saved %d files\n", length(results)))