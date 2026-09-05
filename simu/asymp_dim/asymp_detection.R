# ===============================================
# asymp_detection on Great Lakes -- compute only, no plotting.
# Single job, within-script parallelization over the (beta, r, type) grid.
# Plotting stays local (asymp_detection.R).
# ===============================================
suppressPackageStartupMessages({
  library(matrixStats); library(parallel); library(doSNOW); library(foreach)
})

# ================================================
# Parameters
# ================================================
betas  <- seq(0.1, 0.9, by = 0.1) # Beta values
d      <- 50000                   # Dimension
r_res  <- 20                      # Resolution of the r grid per beta
alpha  <- 0.05                    # Significance level
N      <- as.integer(Sys.getenv("ASYMP_N", "10000"))  # Monte Carlo samples
SMOKE  <- as.integer(Sys.getenv("ASYMP_SMOKE", "0"))  # >0: run only this many grid cells
SEED   <- as.integer(Sys.getenv("ASYMP_SEED", "20260902"))

# chunking to manage memory
chunk_size <- 1000
n_chunks   <- N / chunk_size

outdir <- file.path(getwd(), "asymp_detect_data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
jobid <- Sys.getenv("SLURM_JOB_ID", "local")
tag   <- paste0(jobid, "_", paste(sample(c(letters, 0:9), 6, TRUE), collapse = ""))

cat("N =", N, "| d =", d, "| n_chunks =", n_chunks, "| smoke =", SMOKE, "| tag =", tag, "\n")
flush(stdout())

# ==================================================
# Calculate Null Thresholds (m_alpha and t_alpha)
# ==================================================
cat("Calculating empirical null thresholds...\n"); flush(stdout())
set.seed(SEED)
t0 <- Sys.time()
max_null_dist <- numeric(N)
sum_null_dist <- numeric(N)

for (chk in 1:n_chunks) {
  Y <- matrix(rnorm(chunk_size * d), nrow = chunk_size)
  idx <- ((chk - 1) * chunk_size + 1):(chk * chunk_size)

  max_null_dist[idx] <- rowMaxs(abs(Y))
  P <- 1 / (2 * pnorm(abs(Y), lower.tail = FALSE))
  sum_null_dist[idx] <- rowMeans(P)
  cat(sprintf("  threshold chunk %d/%d (%.1f s)\n", chk, n_chunks,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  flush(stdout())
}

m_alpha <- quantile(max_null_dist, 1 - alpha)
t_alpha <- quantile(sum_null_dist, 1 - alpha)
cat("Thresholds -> Max Test:", m_alpha, "| Pareto Test:", t_alpha, "\n"); flush(stdout())
rm(Y, P, max_null_dist, sum_null_dist); invisible(gc())

# ====================================================
# Parallel & Dynamic Grid Setup
# ====================================================
# detectCores() reports the whole node, not the cgroup allocation -- always
# take the worker count from Slurm.
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA))
if (is.na(slurm_cpus)) slurm_cpus <- parallel::detectCores()
cores_to_use <- max(min(slurm_cpus - 4, 40), 1)

cl <- parallel::makeCluster(cores_to_use, type = "SOCK")
registerDoSNOW(cl)
parallel::clusterSetRNGStream(cl, SEED)   # reproducible, worker-independent streams
cat("Running parallel simulation on", cores_to_use, "cores (Slurm gave", slurm_cpus, ")\n")
flush(stdout())

# Dynamic parameter grid where r limits depend on beta
param_list <- list()
for (b in betas) {
  k <- round(b / 0.1)

  if (k <= 5) {
    r_max <- 1 + (k - 1) * 0.5     # beta <= 0.5: steps of 0.5
  } else {
    r_max <- 2.0 + (k - 5)    # beta > 0.5: starts at 2.0, steps of 1
  }

  rs_beta <- seq(0, r_max, length.out = r_res)

  for (t in c("I", "II")) {
    param_list[[length(param_list) + 1]] <- data.frame(beta = b, r = rs_beta, type = t)
  }
}
param_grid <- do.call(rbind, param_list)

if (SMOKE > 0) {
  # spread the smoke cells across the grid so the timing is representative
  keep <- unique(round(seq(1, nrow(param_grid), length.out = min(SMOKE, nrow(param_grid)))))
  param_grid <- param_grid[keep, , drop = FALSE]
}
n_cells <- nrow(param_grid)
cat("Grid cells:", n_cells, "\n"); flush(stdout())

# progress: explicit cat + flush, because R buffers stdout when redirected to a file
t_start <- Sys.time()
progress <- function(n) {
  el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat(sprintf("[progress] %d/%d cells  elapsed %.1f min  eta %.1f min\n",
              n, n_cells, el, el / n * (n_cells - n)))
  flush(stdout())
}
opts <- list(progress = progress)

# =====================================================================
# Simulation Loop
# =====================================================================
results <- foreach(idx = 1:n_cells, .combine = rbind,
                   .packages = c("matrixStats"),
                   .options.snow = opts) %dopar% {

                     cell_t0 <- Sys.time()
                     b <- param_grid$beta[idx]
                     r_val <- param_grid$r[idx]
                     dist_type <- param_grid$type[idx]

                     d1 <- floor(d^(1 - b))

                     rejections_max <- 0
                     rejections_sum <- 0

                     for (chk in 1:n_chunks) {
                       Z_null <- matrix(rnorm(chunk_size * (d - d1)), nrow = chunk_size)
                       max_null <- rowMaxs(abs(Z_null))
                       sum_null <- rowSums(1 / (2 * pnorm(abs(Z_null), lower.tail = FALSE)))
                       rm(Z_null)

                       Z_sig <- matrix(rnorm(chunk_size * d1), nrow = chunk_size)

                       # Handling r = 0 case explicitly to avoid rcauchy scale=0 (NaN) errors
                       if (r_val == 0) {
                         nu <- matrix(0, nrow = chunk_size, ncol = d1)
                       } else {
                         if (dist_type == "I") {
                           u <- matrix(runif(chunk_size * d1, -0.5, 0.5), nrow = chunk_size)
                           nu <- -r_val * sign(u) * log(1 - 2 * abs(u))
                         } else {
                           scale_cauchy <- (r_val * sqrt(2 * log(d))) / (d^(1 - b))
                           nu <- matrix(rcauchy(chunk_size * d1, location = 0, scale = scale_cauchy), nrow = chunk_size)
                         }
                       }

                       X_sig <- Z_sig + nu
                       rm(Z_sig, nu)

                       max_sig <- rowMaxs(abs(X_sig))
                       sum_sig <- rowSums(1 / (2 * pnorm(abs(X_sig), lower.tail = FALSE)))
                       rm(X_sig)

                       max_val <- pmax(max_null, max_sig)
                       sum_val <- (sum_null + sum_sig) / d

                       rejections_max <- rejections_max + sum(max_val > m_alpha)
                       rejections_sum <- rejections_sum + sum(sum_val > t_alpha)
                     }

                     data.frame(
                       beta = b,
                       r = r_val,
                       type = dist_type,
                       Power_Max = rejections_max / N,
                       Power_Sum = rejections_sum / N,
                       secs = as.numeric(difftime(Sys.time(), cell_t0, units = "secs"))
                     )
                   }

stopCluster(cl)

# =====================================================================
# Save (plotting is done locally)
# =====================================================================
results$N <- N
results$d <- d
results$alpha <- alpha
results$m_alpha <- as.numeric(m_alpha)
results$t_alpha <- as.numeric(t_alpha)
results$seed <- SEED

fout <- file.path(outdir, sprintf("asymp_detection_%s%s.csv", if (SMOKE > 0) "smoke_" else "", tag))
write.csv(results, fout, row.names = FALSE)
cat("Wrote", fout, "with", nrow(results), "rows\n")
cat(sprintf("Total wall time: %.1f min | per-cell secs: median %.1f, max %.1f\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            median(results$secs), max(results$secs)))
flush(stdout())
