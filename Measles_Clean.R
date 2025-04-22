# rm(list = ls())

#' Jonathan P. Smith, PhD, MPH
#' Yale University
#  jonathan.p.smith@yale.edu
#' Last updated: April 22, 2025


#' Step 1: Read in and prepare data from U.S. CDC
#' Data from https://www.cdc.gov/measles/data-research/index.html

filep <- "~jonathansmith/Downloads/"
xx <- read.csv(paste0(filep, "data-table-3.csv"))
xx <- xx[grepl("2000-Present", xx$filter),] # cases before 2000 do no have proportion of zeros
xx$p <- (xx$cases - xx$outbreaks_cases) / xx$cases

#' Step 2 quantify k using proportion of zeros method from Lloyd-Smith et al
#' Supplementary Notes section 2.2.2, https://static-content.springer.com/esm/art%3A10.1038%2Fnature04153/MediaObjects/41586_2005_BFnature04153_MOESM1_ESM.pdf

#' function quantifies k for a given value of R and p:
#'  @param R Basic reproducive number
#'  @param p proportion of individuals who did not transmit any secondary cases
solve_k <- Vectorize(function(R, p, lower = 1e-6, upper = 1e4) {
  f <- function(k)
    (1 + R / k)^(-k) - p
  tryCatch({
    uniroot(f, lower = lower, upper = upper)$root
  }, error = function(e) Inf)
})


#' Step 3 - define function for extinction probabilities
#' @param R basic reproduction number
#' @param k dispersion parameter
extinction_prob <- function(R, k, tol = 1e-8) {
  if (R <= 1) # probability of extinction is 1 if R<1
    (1)
  G <- \(q) (k / (k + R * (1 - q)))^k
  f <- \(q) G(q) - q
  uniroot(f, lower = 0, upper = 1 - tol)$root
}

# Step 4 - generate background data for heatmap and summary estimates
# Grid of R and k
R_vals <- seq(1.01, 20, by = 0.01)
k_vals <- seq(0.01, 1.5, by = 0.01)
ext_matrix <- matrix(NA, nrow = length(k_vals), ncol = length(R_vals))
for (i in seq_along(k_vals)) {
  for (j in seq_along(R_vals)) {
    ext_matrix[i, j] <- extinction_prob(R_vals[j], k_vals[i])
  }
}
zvals <- t(ext_matrix)

# Uncertainty - 95% CI for 2025 case counts
ci_2025 <- 1 - binom.test(xx$outbreaks_cases[xx$year %in% 2025], xx$cases[xx$year %in% 2025])$conf.int  # returns lower and upper bounds

# Uncertainty - median and IQR for 2000 - 2024 
summary_all <- apply(sapply(xx$p[xx$year < 2025], \(x) solve_k(R_vals, x)), 1, summary)

# Step 5: Plot
# Color map and breaks
n_colors <- 100
colormap <- rev(heat.colors(n_colors))
zlim <- range(zvals, na.rm = TRUE)
breaks <- seq(zlim[1], zlim[2], length.out = n_colors + 1)

# for arrows
xpos <- 5.25

# layout and organize for plots
layout(matrix(c(1,2), nrow = 1), widths = c(5,1))
par(mar = c(5, 4, 4, 0.5))

# Plot heatmap base
image(x = R_vals, y = k_vals, z = zvals, col = colormap, breaks = breaks,
      xlab = "Basic Reproduction Number, R", ylab = "Dispersion parameter, k", useRaster = TRUE,
      axes = FALSE)
axis(1, at =  c(min(R_vals), unique(as.integer(R_vals))[-1]),
     label = unique(as.integer(R_vals)))
axis(2, at = unique(round(k_vals, 1)), las = 2) 

# Overlay 2025 Data
lines(x = R_vals, y = solve_k(R_vals, xx$p[xx$year %in% 2025]), 
      lty = 1, lwd = 2)
polygon(x = c(R_vals, rev(R_vals)),
        y = c(solve_k(R_vals, ci_2025[1]), rev(solve_k(R_vals, ci_2025[2]))),
        col = adjustcolor("black", alpha.f = 0.2), border = NA)

# Overlay median and IQR for 2000-2024
lines(x = R_vals, y = summary_all[3,],
      lty = 2, lwd = 2)
polygon(x = c(R_vals, rev(R_vals)),
        y = c(summary_all[2,], rev(summary_all[5,])),
        col = adjustcolor("black", alpha.f = 0.2), border = NA)

# Add in arrows and descriptive text
arrows(x0 = xpos, y0 = 1.05, y1 = 1.25, length = 0.1)
text("Transmission\nis more uniform,\nsuggesting sporatic\nextinction is less likely.", 
     x = xpos+ 0.250, y = 1.15, cex = 0.8,
     pos = 4)
arrows(x0 = xpos, y0 = 0.95, y1 = 0.75, length = 0.1)
text("Transmission\nis more heterogeneous,\nsuggesting sporatic extinction \nis more likely.",
     x = xpos + 0.250, y = 0.85, cex = 0.8,
     pos = 4)

# Add legend
legend("topright", c("2025 Estimates", "Median Estimates (2000-2024)", "Uncertainty"),
       #title = "U.S. Measles Outbreaks",
       col = c(1, 1, adjustcolor("black", alpha.f = 0.2)),
       lty = c(1:2, NA), pch = c(NA, NA, 22), 
       pt.cex = 2, pt.bg = c(NA, NA, adjustcolor("black", alpha.f = 0.2)),
       bty = "n")

# Plot legend
par(mar = c(5, 1, 4, 5))
legend_y <- seq(zlim[1], zlim[2], length.out = n_colors)
legend_z <- matrix(legend_y, nrow = 1)  # 1 row, 100 columns

image(x = c(1, 2), y = legend_y, z = legend_z,
      col = rev(colormap), xlab = "", ylab = "", axes = FALSE)
axis(4, las = 1, at = seq(0, 1, 0.1), labels = rev(seq(0, 1, 0.1)))
mtext("Extinction\nProbability", side = 3)#, line = 2.5)

