# ===============================================
# Plot empirical power curves from the Great Lakes run of
# asymp_detection-cluster.R.  Reads the pulled CSV; no simulation here.
# Panels: distribution (rows) x gamma (columns).
# The simulation calls the sparsity parameter `beta`; it is renamed `gamma`
# for plotting.
# ===============================================
library(lattice)
library(latticeExtra)

datadir <- file.path(getwd(), "asymp_detect_data")
figdir  <- file.path(dirname(getwd()), "figs")

# newest full-grid CSV (smoke files are excluded)
csvs <- list.files(datadir, pattern = "^asymp_detection_[0-9]", full.names = TRUE)
stopifnot(length(csvs) > 0)
infile <- csvs[which.max(file.mtime(csvs))]
cat("Reading", infile, "\n")

results <- read.csv(infile, stringsAsFactors = FALSE)
names(results)[names(results) == "beta"] <- "gamma"

gammas <- sort(unique(results$gamma))

# =====================================================================
# Long format: one row per (gamma, r, type, Test)
# =====================================================================
res_max <- results[, c("gamma", "r", "type", "Power_Max")]
res_max$Test <- "Max test"
names(res_max)[4] <- "Power"

res_sum <- results[, c("gamma", "r", "type", "Power_Sum")]
res_sum$Test <- "Pareto combination test"
names(res_sum)[4] <- "Power"

plot_data <- rbind(res_max, res_sum)

plot_data$Test  <- factor(plot_data$Test, levels = c("Pareto combination test", "Max test"))
plot_data$gammaF <- factor(plot_data$gamma, levels = gammas)
plot_data$distF  <- factor(ifelse(plot_data$type == "I", "Laplace", "Cauchy"),
                           levels = c("Laplace", "Cauchy"))

gamma_expr <- as.expression(lapply(gammas, function(g) bquote(gamma == .(g))))

# =====================================================================
# Single figure faceted on distribution and gamma
# =====================================================================
p <- xyplot(Power ~ r | gammaF + distF,
            data = plot_data,
            groups = Test,
            type = "l",

            # h/v = -1 puts the grid on the axis ticks, so it follows the
            # per-panel free x-scale
            panel = function(...) {
              panel.grid(h = 3, v = 3, lty = 3, col = "grey80")
              panel.abline(h=c(0,1), lty=3, col = "grey80")
              panel.xyplot(...)
            },

            auto.key = list(columns = 2, space = "top", lines = TRUE, points = FALSE),

            # order follows levels(plot_data$Test): Pareto first, then Max
            par.settings = list(
              superpose.line = list(col = c("red", "blue3"), lty = c(1, 2), lwd = 2)
            ),

            # tick.number keeps neighbouring panels' free x-labels from colliding
            scales = list(x = list(relation = "free", tick.number = 3, cex = 0.7),
                          y = list(cex = 0.7),
                          alternating = 1),
            as.table = TRUE,   # Laplace row on top, not bottom
            xlab = "Signal strength r",
            ylab = "Empirical power",
            ylim = c(-0.05, 1.05))

# gamma on the top strips, distribution on the left strips
p <- useOuterStrips(p,
                    strip = strip.custom(factor.levels = gamma_expr),
                    strip.left = strip.custom(horizontal = FALSE))
print(p)
