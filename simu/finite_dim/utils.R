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

# combination tests -----
combine.test <- function(P.vals, w=NULL, method="Pareto") {
  stopifnot(method %in% c("Pareto", "Cauchy", "Cauchy+","Frechet"))
  d <- dim(P.vals)[2]
  n <- dim(P.vals)[1]
  if (is.null(w)) {
    # equal weights
    w <- rep(1/d, d)
  } else {
    stopifnot(all(w > 0))
    w <- w / sum(w)  }
  # uniform scale
  if (method=="Pareto") {
    Y <- 1 / (P.vals)
    Y.combined <- c(Y %*% w)
    pval <- 1 / Y.combined
  } else if (method=="Cauchy") {
    Y <- qcauchy(1-P.vals)
    Y.combined <- c(Y %*% w)
    pval <- 1 - pcauchy(Y.combined)
  } else if (method=="Cauchy+") {
    # absolute value of a Cauchy
    Y <- tan(pi * (1-P.vals) / 2)
    Y.combined <- c(Y %*% w)
    pval <- 1 - 2 / pi * atan(Y.combined)
  } else if (method=="Frechet") {
    P.vals[P.vals == 0] <- 1e-18
    P.vals[P.vals == 1] <- 1 - 1e-18
    # -1 / log(1 - P.vals)
    Y <- -1 / log1p(-P.vals)
    Y.combined <- apply(Y %*% diag(w), 1, max)
    pval <- 1 - exp(-1 / Y.combined)
  }
  return(pval)
}
