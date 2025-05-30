set.seed(20241028)
library(rethinking)
library(plot.matrix)
library(viridis)

L <- 40
data <- expand.grid(X = 1:L, Y = 1:L)
N_points <- L^2

mat <- matrix(0, nrow = L, ncol = L)
center1 <- c(L*0.35, L*0.35)
radius1 <- L*0.05
center2 <- c(L*0.7, L*0.7)
radius2 <- L*0.1

for (i in 1:L) {
  for (j in 1:L) {
    distance1 <- sqrt((i - center1[1])^2 + (j - center1[2])^2)
    distance2 <- sqrt((i - center2[1])^2 + (j - center2[2])^2)
    if (distance1 <= radius1) mat[i, j] <- 1
    if (distance2 <= radius2) mat[i, j] <- 1
  }
}

x <- as.vector(t(mat))
m <- as.matrix(dist(data, method = "euclidean"))

beta <- 5
rho <- sqrt(0.5)
etasq_vals <- c(10, 4, 2, 1, 0.1)
noise_sd <- c(3, 2, 1.5, 1, 0.5)
corr_labels <- c("Random", "Weak", "Moderate", "Medium", "High")

sim_maldi_list <- list()

for(i in 1:5) {
  K <- etasq_vals[i] * exp(-0.5 * ((m / rho)^2)) + diag(1e-9, N_points)
  sim_gp <- MASS::mvrnorm(1, mu = rep(0, N_points), Sigma = K)
  sim_maldi_list[[i]] <- rnorm(N_points, sim_gp + beta * x, sd = noise_sd[i])
}

vec_vessels <- rep(x, 5)
vec_maldis <- unlist(sim_maldi_list)

hilbert_curve <- function(n) {
  if (n == 0) return(matrix(1, 1, 1))
  
  h_prev <- hilbert_curve(n - 1)
  size <- 2^(n-1)
  
  h1 <- t(h_prev[size:1, ])
  h2 <- h_prev + size^2
  h3 <- h_prev + 2 * size^2
  h4 <- t(h_prev)[size:1, size:1] + 3 * size^2
  
  top <- cbind(h1, h2)
  bottom <- cbind(h4, h3)
  rbind(top, bottom)
}

n_hilbert <- ceiling(log2(L))
L_hilbert <- 2^n_hilbert
h_matrix <- hilbert_curve(n_hilbert)

hilbert_order <- numeric(N_points)
idx <- 1
for (i in 1:L) {
  for (j in 1:L) {
    hilbert_order[idx] <- h_matrix[i, j]
    idx <- idx + 1
  }
}

mat_vector <- as.vector(t(mat))
hilbert_vessel <- mat_vector[order(hilbert_order)]

hilbert_vessels <- rep(hilbert_vessel, 5)
hilbert_maldis <- numeric(5 * N_points)

for(k in 1:5) {
  maldi_vector <- sim_maldi_list[[k]]
  hilbert_maldis[((k-1)*N_points + 1):(k*N_points)] <- maldi_vector[order(hilbert_order)]
}

corr_factor <- factor(rep(corr_labels, each = N_points), levels = corr_labels)

lm_vec <- lm(vec_maldis ~ vec_vessels * corr_factor)
lm_hilbert <- lm(hilbert_maldis ~ hilbert_vessels * corr_factor)

beta_vec <- coef(lm_vec)[grep("vec_vessels", names(coef(lm_vec)))]
beta_vec <- c(beta_vec[1], beta_vec[1] + beta_vec[-1])
names(beta_vec) <- corr_labels

beta_hilbert <- coef(lm_hilbert)[grep("hilbert_vessels", names(coef(lm_hilbert)))]
beta_hilbert <- c(beta_hilbert[1], beta_hilbert[1] + beta_hilbert[-1])
names(beta_hilbert) <- corr_labels

se_vec <- summary(lm_vec)$coefficients[grep("vec_vessels", rownames(summary(lm_vec)$coefficients)), "Std. Error"]
se_hilbert <- summary(lm_hilbert)$coefficients[grep("hilbert_vessels", rownames(summary(lm_hilbert)$coefficients)), "Std. Error"]

if (length(se_vec) < 5) se_vec <- rep(se_vec[1], 5)
if (length(se_hilbert) < 5) se_hilbert <- rep(se_hilbert[1], 5)

dev.new(width = 10, height = 12)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

for(i in 1:5) {
  plot(matrix(sim_maldi_list[[i]], ncol = L, byrow = TRUE), col = viridis,
       key = NULL, main = paste(corr_labels[i], "Correlation"),
       xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
}

plot(mat, col = viridis, key = NULL, main = "Vessel Pattern",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)

ylim_range <- range(c(beta_vec - 2*se_vec, beta_hilbert + 2*se_hilbert), na.rm = TRUE)
if (!all(is.finite(ylim_range))) {
  ylim_range <- c(0, 10)
}

dev.new(width = 8, height = 6)
plot(1:5, beta_vec, type = "n", ylim = ylim_range,
     xlab = "Correlation Level", ylab = "Beta Coefficient", xaxt = "n",
     main = "Regression Coefficients: Vector vs Hilbert")
axis(1, at = 1:5, labels = corr_labels)
abline(h = beta, lty = 2, col = "gray50")

points(1:5 - 0.1, beta_vec, pch = 19, col = "blue", cex = 1.5)
segments(1:5 - 0.1, beta_vec - 1.96*se_vec, 1:5 - 0.1, beta_vec + 1.96*se_vec, col = "blue", lwd = 2)

points(1:5 + 0.1, beta_hilbert, pch = 19, col = "red", cex = 1.5)
segments(1:5 + 0.1, beta_hilbert - 1.96*se_hilbert, 1:5 + 0.1, beta_hilbert + 1.96*se_hilbert, col = "red", lwd = 2)

legend("topright", legend = c("Vectorized", "Hilbert", "True Beta"), 
       col = c("blue", "red", "gray50"), pch = c(19, 19, NA), lty = c(NA, NA, 2), cex = 0.9)