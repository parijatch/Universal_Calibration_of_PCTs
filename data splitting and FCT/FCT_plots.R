library(plyr)
library(MASS)
library(data.table)
library(latex2exp)

source("../utils.R")

######Simulate Multivariate-t######
simu.mvt.null <- function(d=4, n=1e5, nu=1, cor.mat=diag(d)){
  
  stopifnot(all(diag(cor.mat)==1))
  # multivariate t with cov=cor.mat
  X <- mvrnorm(n=n, mu=rep(0, d), Sigma=cor.mat)
  denom <- replicate(n, sqrt(rgamma(n=1,shape = nu/2,rate = 1/2)/nu))
  X <- X / denom
  
  pval <- 1 - pt(X, df=nu)
  return(pval)
}

########Generate Splits#########
generate_subsets <- function(d, k, min_size = 2, max_size = NULL) {
  if (is.null(max_size)) {
    max_size <- ceiling(d / 2)
  }
  
  subsets <- list()
  uncovered <- 1:d
  j <- 1
  
  # Fill up to k - 1 subsets or until all indices are covered
  while (length(uncovered) > 0 && j < k) {
    size_j <- sample(min_size:max_size, 1)
    subset_j <- sort(sample(1:d, size_j, replace = FALSE))
    subsets[[j]] <- subset_j
    uncovered <- setdiff(uncovered, subset_j)
    j <- j + 1
  }
  
  # Final subset: must include all remaining uncovered elements
  remaining <- setdiff(1:d, unlist(subsets))
  subset_min <- max(min_size, length(remaining))
  subset_max <- max(subset_min, max_size)
  size_k <- sample(subset_min:subset_max, 1)
  
  extra_pool <- setdiff(1:d, remaining)
  extra_indices <- if (length(extra_pool) >= (size_k - length(remaining))) {
    sample(extra_pool, size_k - length(remaining), replace = FALSE)
  } else {
    sample(extra_pool, size_k - length(remaining), replace = TRUE)
  }
  
  subsets[[k]] <- sort(unique(c(remaining, extra_indices)))
  return(subsets)
}


######Extract sidak corrected p value matrix######
transform_pvalues <- function(P1, subsets) {
  n <- nrow(P1)
  k <- length(subsets)
  P2 <- matrix(NA, nrow = n, ncol = k)
  
  for (j in 1:k) {
    cols <- subsets[[j]]
    min_vals <- apply(P1[, cols, drop = FALSE], 1, min)
    P2[, j] <- 1 - (1 - min_vals)^length(cols)
  }
  
  return(P2)
}

######Calibration Lineplot######
calibration.lineplot.split <- function(nu.vec = c(1, 30),
                                 alpha.vec = 10^seq(log10(0.001), log10(0.1), length.out = 100),
                                 d = 10, n = 1e6,k=5,rho=0.5,
                                 cor.type = "autoreg",
                                 methods = c("Pareto", "Cauchy", "Frechet")) {
  stopifnot(all(methods %in% c("Pareto", "Cauchy", "Frechet")))
  
  method_colors <- c(Pareto = "red", Cauchy = "blue", Frechet = "darkgreen")
  method_labels <- c(Pareto = "PCT", Cauchy = "CCT", Frechet = "FCT")
  
  par(mfrow = c(1, length(nu.vec)), mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))
  
  for (nu in nu.vec) {
    message(sprintf("Running nu = %f", nu))
    cor.mat <- get.cor(d, rho=rho,cor.type = cor.type)
    
    # simulate under null
    
    pval<-simu.mvt.null(d=d,n=n,nu=nu,cor.mat = cor.mat)
    splits<-generate_subsets(d=dim(pval)[2],k=k)
    pval.sidak<-transform_pvalues(pval,splits)
    
    # compute rejection rates
    rejections <- sapply(alpha.vec, function(alpha) {
      result <- numeric(length(methods))
      for (i in seq_along(methods)) {
        method <- methods[i]
        p.combined <- combine.test(pval.sidak, method = method,splits=splits,rej.larger = TRUE)
        result[i] <- mean(p.combined < alpha)
      }
      return(result)
    })
    
    rejections <- t(as.matrix(rejections))
    # rows = alphas, cols = methods
    ratio <- rejections / alpha.vec
    
    lines_to_plot <- if (length(methods) == 1) matrix(ratio, ncol = 1) else ratio
    
    # Plot
    plot(1 / alpha.vec, lines_to_plot[, 1], type = "l", col = method_colors[methods[1]],
         ylim = c(0, max(max(lines_to_plot),1.5)),
         xlab = TeX("$1/\\alpha$"), ylab = "Empirical Rejection Rate / alpha",
         main = TeX(sprintf("ν = %.2f", nu)))
    
    if (length(methods) > 1) {
      for (i in 2:length(methods)) {
        lines(1 / alpha.vec, lines_to_plot[, i], col = method_colors[methods[i]])
      }
    }
    
    
    abline(h = 1, lty = 2)
    legend("bottomright", legend = method_labels[methods],
           col = method_colors[methods], pch = 1)
  }
  
  title_string <- paste(method_labels[methods], collapse = ", ")
  mtext(
    TeX(sprintf("Calibration of %s under Multivariate-$t_\\nu(0,\\Sigma)$ with %s $\\Sigma$",
                title_string, cor.type)),
    outer = TRUE, cex = 1.5, line = 1
  )
  
  par(mfrow = c(1, 1))
}

######Heatmap######
calibration.heatmaps.splits <- function(alpha.vec = c( 0.05, 0.01, 0.005,0.001,0.0005),
                                 nu.vec = c(1, 5, 15, 30),
                                 d = 10, n = 1e6,k=5,rho=0.5,
                                 cor.type = "autoreg") {
  mat.frechet <- matrix(NA, nrow = length(alpha.vec), ncol = length(nu.vec))

  for (i in seq_along(nu.vec)) {
    nu <- nu.vec[i]
    cor.mat <- get.cor(d,rho=rho,cor.type = cor.type)
    
    # simulate null
    pval<-simu.mvt.null(d=d,n=n,nu=nu,cor.mat = cor.mat)
    splits<-generate_subsets(d=dim(pval)[2],k=k)
    pval.sidak<-transform_pvalues(pval,splits)
    
    for (j in seq_along(alpha.vec)) {
      alpha <- alpha.vec[j]
      mat.frechet[j, i] <- mean(combine.test(pval.sidak, method = "Frechet",rej.larger = TRUE,splits = splits) < alpha) / alpha
    }
  }
  
  df.frechet <- melt(mat.frechet)
  
  colnames(df.frechet) <- c("alpha.idx", "nu.idx", "ratio")
  df.frechet$Method <- "Frechet"
  df <- df.frechet
  
  df$alpha <- factor(alpha.vec[df$alpha.idx], levels = rev(alpha.vec))
  df$nu <- factor(nu.vec[df$nu.idx])
  
  gg <- ggplot(df, aes(x = nu, y = alpha, fill = ratio)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 1,limits=c(0.2,1.5),
                         name = TeX("Empirical Rejection / $\\alpha$")) +
    facet_wrap(~Method) +
    labs(title = TeX(sprintf("Calibration Heatmap for FCT with Multivariate-$t_\\nu(0,\\Sigma)$ with %s $\\Sigma$", cor.type)),
         x = TeX("$\\nu$"), y = TeX("$\\alpha$")) +
    theme_minimal(base_size = 13)
  
  
  print(gg)
}
