set.seed(20241028)
library(rethinking)
library(plot.matrix)
library(viridis)
library(MASS)

# PARAMETERS
num_images <- 3  
n <- 20         

# output folder sul mio desktop
out_folder <- "C:/Users/GIUSEPPE/Desktop/Gaussian Process Model/x"
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

# SIMULATE VESSEL IMAGE
simulate_vessel <- function(n, circles) {
  mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (circle in circles) {
        center <- circle$center
        radius <- circle$radius
        if (sqrt((i - center[1])^2 + (j - center[2])^2) <= radius) {
          mat[i, j] <- 1
        }
      }
    }
  }
  return(mat)
}

# VESSEL DEFINITIONS
vessel_list <- list(
  list( list(center = c(7,7), radius = 1),  list(center = c(14,14), radius = 2) ),
  list( list(center = c(5,15), radius = 1.5), list(center = c(15,5), radius = 2.5) ),
  list( list(center = c(10,10), radius = 2), list(center = c(3,3), radius = 1) )
)

if (length(vessel_list) < num_images) {
  for (i in (length(vessel_list) + 1):num_images) {
    circles <- list(
      list(center = c(sample(1:n, 1), sample(1:n, 1)), radius = runif(1, 0.5, 2)),
      list(center = c(sample(1:n, 1), sample(1:n, 1)), radius = runif(1, 0.5, 2))
    )
    vessel_list[[i]] <- circles
  }
}

# SIMULATE VESSEL IMAGES
vessel_mats <- list()
vessel_vectors <- list()
coords_list <- list()

for (i in 1:num_images) {
  circles <- vessel_list[[i]]
  mat <- simulate_vessel(n, circles)
  vessel_mats[[i]] <- mat
  
  pdf(file.path(out_folder, paste0("simulated_vessel_", i, ".pdf")))
  plot(mat, col = viridis, key = NULL, main = paste("Simulated Vessels", i),
       xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
  dev.off()
  
  vec <- as.vector(t(mat))
  vessel_vectors[[i]] <- vec
  
  grid <- expand.grid(X = 1:n, Y = 1:n)
  grid$X <- grid$X + (i - 1) * n  
  coords_list[[i]] <- grid
}

x_combined <- do.call(c, vessel_vectors)
data_combined <- do.call(rbind, coords_list)
colnames(data_combined) <- c("X_Coordinate", "Y_Coordinate")
N_points <- nrow(data_combined)

m <- as.matrix(dist(data_combined, method = "euclidean"))

# PRIOR KERNEL VISUALIZATION
p.etasq <- rexp(N_points, 2)[1:20]
p.rhosq <- rexp(N_points, 0.5)[1:20]
pdf(file.path(out_folder, "kernel_prior.pdf"))
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 1),
     xlab = "Pixel distance", ylab = "Covariance", main = "Prior Simulation")
for (i in 1:20) {
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
}
dev.off()

# GP SIMULATION
beta <- 5
etasq <- 2
rho <- sqrt(0.5)

K <- etasq * exp(-0.5 * ((m / rho)^2)) + diag(1e-9, N_points)
sim_gp <- MASS::mvrnorm(1, mu = rep(0, N_points), Sigma = K)
sim_y <- rnorm(N_points, mean = sim_gp + beta * x_combined, sd = 1)

# SIMULATED MALDI IMAGES
pdf(file.path(out_folder, "simulated_maldi_combined.pdf"))
par(mfrow = c(1, num_images))
for (i in 1:num_images) {
  im_data <- sim_y[((i - 1) * n * n + 1):(i * n * n)]
  image(matrix(im_data, nrow = n, ncol = n, byrow = TRUE),
        col = viridis, main = paste("Simulated MALDI - Image", i),
        xlab = "", ylab = "", axes = FALSE)
}
dev.off()

# MODEL FITTING
dat_list <- list(
  y = sim_y,
  x = x_combined,
  S = 1:N_points,
  Dmat = m,
  N = N_points
)

GP <- ulam(
  alist(
    y ~ multi_normal(mu, K),
    mu <- a + b * x,
    matrix[N, N]:K <- cov_GPL2(Dmat, etasq, rho, 0.01),
    a ~ normal(0, 1),
    b ~ normal(0, 0.5),
    etasq ~ dexp(2),
    rho ~ dexp(0.5)
  ),
  data = dat_list, chains = 4, cores = 4, iter = 1000
)

precis(GP)
post <- extract.samples(GP)

# KERNEL PLOT
pdf(file.path(out_folder, "kernel_prior_actual_post.pdf"))
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 2),
     xlab = "Pixel distance", ylab = "Covariance", 
     main = "Prior, Actual, and Estimated Kernel")
for (i in 1:20) {
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
}
curve(etasq * exp(-0.5 * rho * x^2), add = TRUE, lwd = 4)
for (i in 1:20) {
  curve(post$etasq[i] * exp(-post$rho[i] * x), add = TRUE,
        col = col.alpha(4, 0.3), lwd = 6)
}
legend("topright", lwd = 4, col = c(2, 1, 4),
       legend = c("Prior", "Actual", "Estimated"))
dev.off()
