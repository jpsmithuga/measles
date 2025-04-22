# rm(list = ls())
R <- 5:25
#p2025 <- (712 - 660) / 712
tcases_17APR2025 <- 800
ccases_17APR2025 <- 750

p2025 <- (800 - 750) / 800
p2024 <- (285 - 198) / 285
p2023 <- (59 - 29) / 59

solve_k <- Vectorize(function(R, p, lower = 1e-6, upper = 1e4) {
  f <- function(k)
    (1 + R / k)^(-k) - p
  uniroot(f, lower = lower, upper = upper)$root
})

ci_p <- 1 - binom.test(ccases_17APR2025, tcases_17APR2025)$conf.int  # returns lower and upper bounds
ci_p_2024 <- 1 - binom.test(198, 285)$conf.int  # returns lower and upper bounds
ci_p_2023 <- 1 - binom.test(29, 59)$conf.int  # returns lower and upper bounds

# ci_matrix <- sapply(R, \(x) solve_k(x, sort(c(p2025, ci_p))))
# ci_matrix2024 <- sapply(R, \(x) solve_k(x, sort(c(p2024, ci_p_2024))))

# Extinction probability function
extinction_prob <- function(R, k, tol = 1e-8) {
  if (R <= 1) return(1)
  G <- function(q) (k / (k + R * (1 - q)))^k
  f <- function(q) G(q) - q
  uniroot(f, lower = 0, upper = 1 - tol)$root
}

# Grid of R and k
R_vals <- seq(5, 20, by = 0.1)
k_vals <- seq(0.01, 1.5, by = 0.01)
ext_matrix <- matrix(NA, nrow = length(k_vals), ncol = length(R_vals))

# Fill in matrix
for (i in seq_along(k_vals)) {
  for (j in seq_along(R_vals)) {
    ext_matrix[i, j] <- extinction_prob(R_vals[j], k_vals[i])
  }
}

# Transpose for correct orientation
zvals <- t(ext_matrix)

# Color map
n_colors <- 100
colormap <- rev(heat.colors(n_colors))
zlim <- range(zvals, na.rm = TRUE)
breaks <- seq(zlim[1], zlim[2], length.out = n_colors + 1)

# Layout: heatmap + legend
layout(matrix(c(1,2), nrow = 1), widths = c(5,1))

# Plot heatmap
par(mar = c(5, 4, 4, 1))
image(x = R_vals, y = k_vals, z = zvals, col = colormap, breaks = breaks,
      xlab = "Basic Reproduction Number, R", ylab = "Dispersion parameter, k", useRaster = TRUE,
      axes = FALSE)
axis(1, at = unique(as.integer(R_vals)))
axis(2, at = unique(round(k_vals, 1)), las = 2) 

# 2025 
lines(x = R, y = solve_k(R, p2025))
polygon(x = c(R, rev(R)),
        y = c(solve_k(R, ci_p[1]), rev(solve_k(R, ci_p[2]))),
        col = adjustcolor("black", alpha.f = 0.2), border = NA)
# 2024
lines(x = R, y = solve_k(R, p2024), lty = 2)
polygon(x = c(R, rev(R)),
        y = c(solve_k(R, ci_p_2024[1]), rev(solve_k(R, ci_p_2024[2]))),
        col = adjustcolor("black", alpha.f = 0.2), border = NA)
# 2023
lines(x = R, y = solve_k(R, p2023), lty = 3)
polygon(x = c(R, rev(R)),
        y = c(solve_k(R, ci_p_2023[1]), rev(solve_k(R, ci_p_2023[2]))),
        col = adjustcolor("black", alpha.f = 0.2), border = NA)

#Blumberg
# 0.27, 95% CI: 0.18, 0.41
points(x = 4.95, y = 0.27, pch = 23, bg = "darkgrey", xpd = TRUE, cex = 1.5)
# Worden
# 0.40
points(x = 4.95, y = 0.4, pch = 24, bg = "darkgrey", xpd = TRUE, cex = 1.5)

arrows(x0 = min(R)+0.25, y0 = 1, y1 = 1.3, length = 0.1)
text("Transmission\nis more uniform,\nsuggesting sporatic extinction is\nless likely.", x = min(R)+0.50, y = (1+1.3)/2, cex = 0.7,
     pos = 4)

arrows(x0 = min(R)+0.25, y0 = 0.9, y1 = 0.6, length = 0.1)
text("Transmission\nis more heterogeneous,\nsuggesting sporatic extinction \nis more likely.",
     x = min(R)+0.50, y = (0.9+0.6)/2, cex = 0.7,
     pos = 4)

legend("topright", c("2025", "2024", "2023", "95% CI", 
                     "Blumberg, et al (2013)", "Worden, et al (2020)"),
       title = "U.S. Measeles Outbreaks",
       col = c(1,1,1,adjustcolor("black", alpha.f = 0.2), 1, 1),
       lty = c(1:3, NA, NA, NA), pch = c(NA, NA, NA, 22, 23, 24), 
       pt.cex = 1.5, pt.bg = c(NA, NA, NA, 
                             adjustcolor("black", alpha.f = 0.2),
                             "darkgrey", "darkgrey"),
       bty = "n")
# 
# lines(x = R, y = solve_k(R, ci_p_2024[1]), lty = 3)
# lines(x = R, y = solve_k(R, ci_p_2024[2]), lty = 3)

#lines(x = R, y = solve_k(R, p2024))
# polygon(x = c(R, rev(R)),
#         y = c(ci_matrix[1,], rev(ci_matrix[3,])),
#         col = adjustcolor("blue", alpha.f = 0.7), border = NA)
# 
# polygon(x = c(R, rev(R)),
#         y = c(ci_matrix2024[1,], rev(ci_matrix2024[3,])),
#         col = adjustcolor("grey", alpha.f = 0.7), border = NA)


# Plot legend
par(mar = c(5, 1, 4, 5))
legend_y <- seq(zlim[1], zlim[2], length.out = n_colors)
legend_z <- matrix(legend_y, nrow = 1)  # 1 row, 100 columns

image(x = c(1, 2), y = legend_y, z = legend_z,
      col = rev(colormap), xlab = "", ylab = "", axes = FALSE)
axis(4, las = 1, at = seq(0, 1, 0.1), labels = rev(seq(0, 1, 0.1)))
mtext("Extinction\nProbability", side = 3)#, line = 2.5)
