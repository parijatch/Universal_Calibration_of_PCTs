Normal.scale <- function(y, from.uniform=F){
  #
  # Convention: extremes are transformed into extremes
  #
  d = dim(y)[2];
  n = dim(y)[1];
  if (from.uniform==TRUE){
    return(qnorm(y))}
  else{
    x = c();
    for (i in c(1:d)){
      x = cbind(x,qnorm(rank(y[,i],)/(n+1)))
    }
    return(x)
  }
}

Pareto.scale <- function(y,from.uniform=F){
  #
  # Convention: extremes are transformed into extremes
  #
  d = dim(y)[2];
  n = dim(y)[1];
  if (from.uniform==TRUE){
    return(1/(1-y))}
  else{
    x = c();
    for (i in c(1:d)){
      x = cbind(x,(n+1)/rank(-y[,i],))
    }
    return(x)}
}

qfrechet <- function(p, alpha = 1) {
  p <- ifelse(p == 0, 1e-7, ifelse(p == 1, 1 - 1e-7, p))
  (-log(p))^(-1 / alpha)
}

compute_c <- function(n, subsets, w) {
  k <- length(subsets)
  stopifnot(length(w) == k)
  
  c_val <- 0
  for (j in 1:n) {
    max_val <- 0
    for (i in 1:k) {
      if (j %in% subsets[[i]]) {
        val <- w[i] / length(subsets[[i]])
        if (val > max_val) {
          max_val <- val
        }
      }
    }
    c_val <- c_val + max_val
  }
  
  return(c_val)
}

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

combine.test <- function(X, w=NULL, from.uniform=TRUE, rej.larger=TRUE,splits=NULL,
                         method="Pareto") {
  stopifnot(method %in% c("Pareto", "Cauchy", "Cauchy+","Frechet"))
  d <- dim(X)[2]
  n <- dim(X)[1]
  if (is.null(w)) {
    w <- rep(1/d, d)
  } else {
    w <- w / sum(w)  }
  # uniform scale
  if (from.uniform) {
    stopifnot(min(X)>=0 & max(X)<=1)
    if (rej.larger) {
      U <- X
    } else {
      U <- 1 - X
    }
  } else {
    if (rej.larger) {
      rank.mat <- apply(X, 2, rank)
    } else {
      rank.mat <- apply(-X, 2, rank)
    }
    U <- (rank.mat - 1/2) / n
  }
  if (method=="Pareto") {
    Y <- 1 / (U)
    Y.combined <- c(Y %*% w)
    pval <- 1 / Y.combined
  } else if (method=="Cauchy") {
    Y <- qcauchy(1-U)
    Y.combined <- c(Y %*% w)
    pval <- 1 - pcauchy(Y.combined)
  } else if (method=="Cauchy+") {
    # absolute value of a Cauchy
    Y <- tan(pi * (1-U) / 2)
    Y.combined <- c(Y %*% w)
    pval <- 1 - 2 / pi * atan(Y.combined)
  } else if (method=="Frechet") {
    stopifnot(splits!=NULL)
    Y<-qfrechet(1-U)
    scaled_columns <- t(t(Y) * w)  
    Y.combined <- apply(scaled_columns, 1, max)
    c<-compute_c(dim(Y)[1],splits,w)
    Y.combined.scaled<-Y.combined/c
    pval<-1/Y.combined.scaled
  }
  return(pval)
}
