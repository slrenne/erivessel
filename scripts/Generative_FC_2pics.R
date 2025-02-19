set.seed(20241028)
library(rethinking)
library(plot.matrix)
library(viridis)
library(MASS)

# output folder (miodesktop)
out_folder <- "C:/Users/GIUSEPPE/Desktop/Gaussian Process Model/x"
if(!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

# function to simulate a vessel mask matrix
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

# Simulate first vessel image
circles1 <- list(
  list(center = c(7, 7), radius = 1),
  list(center = c(14, 14), radius = 2)
)
mat1 <- test_function(n, circles1)
pdf(file.path(out_folder, "simulated_vessel1.pdf"))
plot(mat1, col = viridis, key = NULL, main = "Simulated vessels 1",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
dev.off()

# Simulate second vessel image
circles2 <- list(
  list(center = c(5, 15), radius = 1.5),
  list(center = c(15, 5), radius = 2.5)
)
mat2 <- test_function(n, circles2)
pdf(file.path(out_folder, "simulated_vessel2.pdf"))
plot(mat2, col = viridis, key = NULL, main = "Simulated vessels 2",
     xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
dev.off()

# Transform matrices into vectors
x1 <- as.vector(t(mat1))
x2 <- as.vector(t(mat2))
x_combined <- c(x1, x2)

# Create coordinate grids
grid1 <- expand.grid(X = 1:n, Y = 1:n)
grid2 <- expand.grid(X = (1:n) + n, Y = 1:n)
data_combined <- rbind(grid1, grid2)
colnames(data_combined) <- c("X_Coordinate", "Y_Coordinate")
N_points <- nrow(data_combined)

# Compute distance matrix
m <- as.matrix(dist(data_combined, method = "euclidean"))
m1 <- as.matrix(dist(mat1, method = "euclidean"))
m2 <- as.matrix(dist(mat2, method = "euclidean"))

# Simulate prior kernel curves
p.etasq <- rexp(N_points, 2)[1:20]
p.rhosq <- rexp(N_points, 0.5)[1:20]
pdf(file.path(out_folder, "kernel_prior.pdf"))
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 1),
     xlab = "pixel distance", ylab = "covariance", main = "Prior simulation")
for(i in 1:20)
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
dev.off()

# Set GP parameters
beta <- 5
etasq <- 2
rho <- sqrt(0.5)

# Generate covariance matrix
K <- etasq * exp(-0.5 * ((m / rho)^2)) + diag(1e-9, N_points)

# Sample latent GP values
sim_gp <- MASS::mvrnorm(1, mu = rep(0, N_points), Sigma = K)

# Generate observed values
sim_y <- rnorm(N_points, mean = sim_gp + beta * x_combined, sd = 1)
sim_y1 <- rnorm(400, sim_gp + beta * x1)
sim_y2 <- rnorm(400, sim_gp + beta * x2)

# Plot simulated MALDI images
pdf(file.path(out_folder, "simulated_maldi_combined.pdf"))
par(mfrow = c(1, 2))
plot(matrix(sim_y[1:(n * n)], nrow = n, ncol = n, byrow = TRUE),
      col = viridis, main = "Simulated MALDI - Image 1",
      xlab = "", ylab = "", axes = FALSE)
plot(matrix(sim_y[(n * n + 1):(2 * n * n)], nrow = n, ncol = n, byrow = TRUE),
      col = viridis, main = "Simulated MALDI - Image 2",
      xlab = "", ylab = "", axes = FALSE)
dev.off()

# Fit Gaussian Process model 2 image separate 

# Prepare data 
dat_list <- list(
  y1 = sim_y1,
  y2 = sim_y2,
  x1 = x1,
  x2 = x2
  # Dmat1 = m1,
  # Dmat2 = m2,
  # N = 20 # side of the matrix
)


GP2 <- ulam(
  alist(
    # GP image 1
    y1 ~ normal(mu1, sigma), 
    mu1 <- a + b * x1,
    #matrix[N, N]:K1 <- cov_GPL2(Dmat1, etasq, rho, 0.01), 
    # GP image 2
    y2 ~ normal(mu2, sigma), 
    mu2 <- a + b * x2,
    #matrix[N, N]:K2 <- cov_GPL2(Dmat2, etasq, rho, 0.01),
    a ~ normal(0, 1),
    b ~ normal(0, 0.5),
    sigma ~ dexp(1)
    #etasq ~ dexp(2),
    #rho ~ dexp(0.5)
  ),   data = dat_list, chains = 4, cores = 4 )

GP2 <- ulam(
  alist(
    # GP image 1
    y1 ~ normal(mu, 0), 
    mu <- a + b * x1,
    #matrix[N, N]:K1 <- cov_GPL2(Dmat1, etasq, rho, 0.01), 
    # GP image 2
    y2 ~ normal(mu, 0), 
    mu <- a + b * x2,
    #matrix[N, N]:K2 <- cov_GPL2(Dmat2, etasq, rho, 0.01),
    a ~ normal(0, 1),
    b ~ normal(0, 0.5),
    #etasq ~ dexp(2),
    #rho ~ dexp(0.5)
  ),   data = dat_list, chains = 1, cores = 1, iter = 100 )




# Prepare data 
dat_list <- list(
  y = sim_y,
  x = x_combined,
  S = 1:N_points,
  Dmat = m,
  N = N_points
)

# Fit Gaussian Process model
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
  data = dat_list, chains = 4, cores = 4, iter = 1000, threads = 2
)

precis(GP)
post <- extract.samples(GP)

# Plot kernel comparison
pdf(file.path(out_folder, "kernel_prior_actual_post.pdf"))
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 2),
     xlab = "pixel distance", ylab = "covariance", main = "Prior, Actual, and Estimated Kernel")
for(i in 1:20)
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
curve(etasq * exp(-0.5 * rho * x^2), add = TRUE, lwd = 4)
for(i in 1:20) {
  curve(post$etasq[i] * exp(-post$rho[i] * x), add = TRUE,
        col = col.alpha(4, 0.3), lwd = 6)
}
legend("topright", lwd = 4, col = c(2, 1, 4),
       legend = c("Prior", "Actual", "Estimated"))
dev.off()
