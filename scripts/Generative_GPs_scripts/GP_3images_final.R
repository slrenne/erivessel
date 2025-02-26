set.seed(20241028)
library(rethinking)
library(plot.matrix)
library(viridis)
library(MASS)

# A helper function to simulate a vessel mask matrix given a list of circles
test_function <- function(n, circles) {
  mat <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      for(circle in circles){
        center <- circle$center
        radius <- circle$radius
        if ( sqrt((i - center[1])^2 + (j - center[2])^2) <= radius ){
          mat[i, j] <- 1
        }
      }
    }
  }
  return(mat)
}

# Set matrix dimensions
n <- 20   # each image will be n x n pixels

# --- Simulate vessel masks for 3 immagini ---
# Image 1: similar to Code 1/2
circles1 <- list(
  list(center = c(7, 7), radius = 1),
  list(center = c(14, 14), radius = 2)
)
mat1 <- test_function(n, circles1)
plot(mat1, col = viridis, key = NULL, main = "Simulated Vessels 1",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)

# Also plot the vector form for image 1 
x1 <- as.vector(t(mat1))
plot(x1, main = "Vector Form of Vessel Mask 1", xlab = "Pixel index", ylab = "Value")

# Image 2 part
circles2 <- list(
  list(center = c(5, 15), radius = 1.5),
  list(center = c(15, 5), radius = 2.5)
)
mat2 <- test_function(n, circles2)
plot(mat2, col = viridis, key = NULL, main = "Simulated Vessels 2",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
x2 <- as.vector(t(mat2))
plot(x2, main = "Vector Form of Vessel Mask 2", xlab = "Pixel index", ylab = "Value")

# Image 3 part
circles3 <- list(
  list(center = c(10, 10), radius = 1.5),
  list(center = c(5, 5), radius = 1)
)
mat3 <- test_function(n, circles3)
plot(mat3, col = viridis, key = NULL, main = "Simulated Vessels 3",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
x3 <- as.vector(t(mat3))
plot(x3, main = "Vector Form of Vessel Mask 3", xlab = "Pixel index", ylab = "Value")

# --- Create coordinate grids for each image ---
# We “shift” the grids so that images do not overlap in the latent process
grid1 <- expand.grid(X = 1:n,       Y = 1:n)
grid2 <- expand.grid(X = (1:n) + n, Y = 1:n)   # shift right by n
grid3 <- expand.grid(X = (1:n) + 2*n, Y = 1:n)   # shift right by 2*n

# Combine the grids
data_combined <- rbind(grid1, grid2, grid3)
N_points <- nrow(data_combined)  # total points (3 * n*n)

# --- Compute the full distance matrix for the combined grid ---
m <- as.matrix(dist(data_combined, method = "euclidean"))

# --- Plot simulated prior kernel curves ---
# (Using only 20 curves as in Code 2)
p.etasq <- rexp(N_points, 2)[1:20]
p.rhosq <- rexp(N_points, 0.5)[1:20]
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 1),
     xlab = "Pixel distance", ylab = "Covariance", main = "Prior Kernel Simulation")
for(i in 1:20){
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6,
        col = col.alpha(2, 0.5))
}

# --- Set GP parameters
beta  <- 5
etasq <- 2
rho   <- sqrt(0.5)

# Generate the covariance matrix using a sqrd exp kernel
K <- etasq * exp(-0.5 * ((m / rho)^2)) + diag(1e-9, N_points)

# Sample 
sim_gp <- MASS::mvrnorm(1, mu = rep(0, N_points), Sigma = K)

# --- Generate MALDI values ---
# Combine the vessel mask vectors (order must match the grid order)
x_combined <- c(x1, x2, x3)
sim_y <- rnorm(N_points, mean = sim_gp + beta * x_combined, sd = 1)

# Split the combined observed data back into three images
sim_y1 <- sim_y[1:(n*n)]
sim_y2 <- sim_y[(n*n + 1):(2*n*n)]
sim_y3 <- sim_y[(2*n*n + 1):(3*n*n)]

# --- Plot the simulated MALDI images ---
par(mfrow = c(1, 3))  # set up a 1x3 plotting layout
image1 <- matrix(sim_y1, nrow = n, ncol = n, byrow = TRUE)
image2 <- matrix(sim_y2, nrow = n, ncol = n, byrow = TRUE)
image3 <- matrix(sim_y3, nrow = n, ncol = n, byrow = TRUE)
plot(image1, col = viridis, key = NULL, main = "Simulated MALDI - Image 1",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
plot(image2, col = viridis, key = NULL, main = "Simulated MALDI - Image 2",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
plot(image3, col = viridis, key = NULL, main = "Simulated MALDI - Image 3",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
par(mfrow = c(1, 1))  # reset layout to one plot per page

# --- Prepare data for fitting three separate GPs with shared parameters ---
dat_list <- list(
  y1    = sim_y1,
  y2    = sim_y2,
  y3    = sim_y3,
  x1    = x1,
  x2    = x2,
  x3    = x3,
  Dmat1 = as.matrix(dist(grid1, method = "euclidean")),
  Dmat2 = as.matrix(dist(grid2, method = "euclidean")),
  Dmat3 = as.matrix(dist(grid3, method = "euclidean")),
  n_sq  = n * n   # number of pixels per image
)

# --- Fit the 3-GP ---
# common regression (a + b*x1,2,3)
GP3 <- ulam(
  alist(
    # Image 1
    y1 ~ multi_normal( mu1 , K1 ),
    mu1 <- a + b * x1,
    matrix[n_sq, n_sq]:K1 <- cov_GPL2( Dmat1, etasq, rho, 0.01 ),
    
    # Image 2
    y2 ~ multi_normal( mu2 , K2 ),
    mu2 <- a + b * x2,
    matrix[n_sq, n_sq]:K2 <- cov_GPL2( Dmat2, etasq, rho, 0.01 ),
    
    # Image 3
    y3 ~ multi_normal( mu3 , K3 ),
    mu3 <- a + b * x3,
    matrix[n_sq, n_sq]:K3 <- cov_GPL2( Dmat3, etasq, rho, 0.01 ),
    
    # Shared parameters 
    a      ~ normal(0, 1),
    b      ~ normal(0, 0.5),
    etasq  ~ dexp(2),
    rho    ~ dexp(0.5)
  ),
  data    = dat_list,
  chains  = 4,
  cores   = 4,
  iter    = 300 #i put 300 to make it compile (around 9 min for 3 images) 
)

precis(GP3)
post <- extract.samples(GP3)

# --- Plot the Prior, Actual, and Estimated Kernel Curves ---
# Here, "Actual" kernel is based on the true values used in simulation.
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 2),
     xlab = "Pixel distance", ylab = "Covariance",
     main = "Prior, Actual, and Estimated Kernel")
# Prior curves (first 20 simulations)
for(i in 1:20){
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6,
        col = col.alpha(2, 0.5))
}
# Actual kernel curve
curve(etasq * exp(-0.5 * (x / rho)^2), add = TRUE, lwd = 4, col = "black")
# Estimated kernel curves (plot a subset from the posterior samples)
n_est <- min(20, length(post$etasq))
for(i in 1:n_est){
  curve(post$etasq[i] * exp(-0.5 * (x / post$rho[i])^2), add = TRUE, 
        col = col.alpha(4, 0.3), lwd = 6)
}
legend("topright", lwd = 4, col = c(2, "black", 4),
       legend = c("Prior", "Actual", "Estimated"))
