set.seed(20241028)
library(rethinking)
library(plot.matrix)
library(viridis)
library(MASS)

# function to simulate a vessel mask matrix (reusing from code 2)
test_function <- function(n, circles) {
  mat <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      for(circle in circles) {
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

# Set matrix dimensions
n <- 20

# Simulate three vessel images with different patterns
circles1 <- list(
  list(center = c(7, 7), radius = 1),
  list(center = c(14, 14), radius = 2)
)

circles2 <- list(
  list(center = c(5, 15), radius = 1.5),
  list(center = c(15, 5), radius = 2.5)
)

circles3 <- list(
  list(center = c(10, 10), radius = 2),
  list(center = c(4, 16), radius = 1.8)
)

# Generate matrices
mat1 <- test_function(n, circles1)
mat2 <- test_function(n, circles2)
mat3 <- test_function(n, circles3)

# Plot vessel images
par(mfrow = c(1, 3))
plot(mat1, col = viridis, key = NULL, main = "Simulated vessels 1",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
plot(mat2, col = viridis, key = NULL, main = "Simulated vessels 2",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
plot(mat3, col = viridis, key = NULL, main = "Simulated vessels 3",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)

# Transform matrices into vectors
x1 <- as.vector(t(mat1))
x2 <- as.vector(t(mat2))
x3 <- as.vector(t(mat3))
x_combined <- c(x1, x2, x3)

# Plot vectors
par(mfrow = c(1, 3))
plot(x1, main = "Vector form - Image 1", type = "l")
plot(x2, main = "Vector form - Image 2", type = "l")
plot(x3, main = "Vector form - Image 3", type = "l")

# Create coordinate grids for each image
grid1 <- expand.grid(X = 1:n, Y = 1:n)
grid2 <- expand.grid(X = (1:n) + n, Y = 1:n)
grid3 <- expand.grid(X = (1:n) + 2*n, Y = 1:n)
data_combined <- rbind(grid1, grid2, grid3)
colnames(data_combined) <- c("X_Coordinate", "Y_Coordinate")
N_points <- nrow(data_combined)

# Compute distance matrices for each image
m <- as.matrix(dist(data_combined, method = "euclidean"))
m1 <- as.matrix(dist(grid1, method = "euclidean"))
m2 <- as.matrix(dist(grid2, method = "euclidean"))
m3 <- as.matrix(dist(grid3, method = "euclidean"))

# Simulate prior kernel curves (following code 1)
p.etasq <- rexp(N_points, 2)[1:20]
p.rhosq <- rexp(N_points, 0.5)[1:20]


# Plot prior kernel curves
par(mfrow = c(1, 1))
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 1),
     xlab = "pixel distance", ylab = "covariance", main = "Prior simulation")
for(i in 1:20)
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6, col = col.alpha(2, 0.5))

# Set GP parameters
beta <- 5
etasq <- 2
rho <- sqrt(0.5)

# Generate covariance matrix
K <- etasq * exp(-0.5 * ((m / rho)^2)) + diag(1e-9, N_points)

# Sample latent GP values
sim_gp <- MASS::mvrnorm(1, mu = rep(0, N_points), Sigma = K)

# Generate observed values for each image
sim_y1 <- rnorm(n*n, mean = sim_gp[1:(n*n)] + beta * x1, sd = 1)
sim_y2 <- rnorm(n*n, mean = sim_gp[(n*n + 1):(2*n*n)] + beta * x2, sd = 1)
sim_y3 <- rnorm(n*n, mean = sim_gp[(2*n*n + 1):(3*n*n)] + beta * x3, sd = 1)

# Plot simulated MALDI images
par(mfrow = c(1, 3))
plot(matrix(sim_y1, nrow = n, ncol = n, byrow = TRUE),
     col = viridis, main = "Simulated MALDI - Image 1",
     xlab = "", ylab = "", axes = FALSE)
plot(matrix(sim_y2, nrow = n, ncol = n, byrow = TRUE),
     col = viridis, main = "Simulated MALDI - Image 2",
     xlab = "", ylab = "", axes = FALSE)
plot(matrix(sim_y3, nrow = n, ncol = n, byrow = TRUE),
     col = viridis, main = "Simulated MALDI - Image 3",
     xlab = "", ylab = "", axes = FALSE)

# Prepare data for fitting
dat_list <- list(
  y1 = sim_y1,
  y2 = sim_y2,
  y3 = sim_y3,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  Dmat1 = m1,
  Dmat2 = m2,
  Dmat3 = m3,
  N = n*n  # number of points per image
)

# Fit three GPs with shared parameters
GP3 <- ulam(
  alist(
    # GP for image 1
    y1 ~ normal(mu1, sigma),
    mu1 <- a + b * x1,
    matrix[N,N]:K1 <- cov_GPL2(Dmat1, etasq, rho, 0.01),
    
    # GP for image 2
    y2 ~ normal(mu2, sigma),
    mu2 <- a + b * x2,
    matrix[N,N]:K2 <- cov_GPL2(Dmat2, etasq, rho, 0.01),
    
    # GP for image 3
    y3 ~ normal(mu3, sigma),
    mu3 <- a + b * x3,
    matrix[N,N]:K3 <- cov_GPL2(Dmat3, etasq, rho, 0.01),
    
    # Shared parameters
    a ~ normal(0, 1),
    b ~ normal(0, 0.5),
    sigma ~ dexp(1),
    etasq ~ dexp(2),
    rho ~ dexp(0.5)
  ), data = dat_list, chains = 4, cores = 4, iter = 300)

# Extract and analyze results
print(precis(GP3))
post <- extract.samples(GP3)

# Plot kernel comparison
par(mfrow = c(1, 1))
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 2),
     xlab = "pixel distance", ylab = "covariance",
     main = "Prior, Actual, and Estimated Kernel")
# Plot priors
for(i in 1:20)
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2),
        add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
# Plot actual kernel
curve(etasq * exp(-0.5 * rho * x^2), add = TRUE, lwd = 4)
# Plot estimated kernels
for(i in 1:20) {
  curve(post$etasq[i] * exp(-post$rho[i] * x),
        add = TRUE, col = col.alpha(4, 0.3), lwd = 6)
}
legend("topright", lwd = 4, col = c(2,1,4),
       legend = c("Prior", "Actual", "Estimated"))
