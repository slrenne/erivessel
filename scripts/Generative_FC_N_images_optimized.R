# Install required packages if not already installed
if (!require("rethinking")) install.packages("rethinking")
if (!require("plot.matrix")) install.packages("plot.matrix")
if (!require("viridis")) install.packages("viridis")
if (!require("MASS")) install.packages("MASS")
if (!require("Rcpp")) install.packages("Rcpp")
if (!require("RcppArmadillo")) install.packages("RcppArmadillo")

# Load libraries
library(rethinking)
library(plot.matrix)
library(viridis)
library(MASS)
library(Rcpp)
library(RcppArmadillo)

# Random seed
set.seed(20241028)

# C++ function for distance matrix
cppFunction(depends = "RcppArmadillo",
'arma::mat computeDistanceMatrix(const arma::mat& coordinates) {
    int n = coordinates.n_rows;
    arma::mat dist_mat(n, n);
    for(int i = 0; i < n; i++) {
        for(int j = i; j < n; j++) {
            double d = sqrt(sum(pow(coordinates.row(i) - coordinates.row(j), 2)));
            dist_mat(i,j) = d;
            dist_mat(j,i) = d;
        }
    }
    return dist_mat;
}')

# output folder sul mio desktop
out_folder <- "C:/Users/GIUSEPPE/Desktop/Gaussian Process Model/x"
if(!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

# Vessel mask function
test_function <- function(n, circles) {
  x_coords <- matrix(rep(1:n, n), nrow = n)
  y_coords <- t(x_coords)
  mat <- matrix(0, nrow = n, ncol = n)
  for(circle in circles) {
    center <- circle$center
    radius <- circle$radius
    distances <- sqrt((x_coords - center[1])^2 + (y_coords - center[2])^2)
    mat <- mat + (distances <= radius)
  }
  return((mat > 0) * 1)
}

# Generate random circles
generate_random_circles <- function(n_circles) {
  circles <- list()
  for(i in 1:n_circles) {
    circles[[i]] <- list(
      center = c(sample(5:15, 1), sample(5:15, 1)),
      radius = runif(1, 1, 2.5)
    )
  }
  return(circles)
}

# Parameters
n <- 20  
n_images <- 5  
n_circles_per_image <- 2  

# Generate vessel images
matrices <- list()
x_all <- c()

for(img in 1:n_images) {
  circles <- generate_random_circles(n_circles_per_image)
  mat <- test_function(n, circles)
  matrices[[img]] <- mat
  
  pdf(file.path(out_folder, sprintf("simulated_vessel%d.pdf", img)))
  plot(mat, col = viridis, key = NULL, 
       main = sprintf("Simulated vessels %d", img),
       xlab = "", ylab = "", axis.col = NULL, axis.row = NULL)
  dev.off()
  
  x_all <- c(x_all, as.vector(t(mat)))
}

# Coordinate grid
grid_list <- list()
for(img in 1:n_images) {
  grid_list[[img]] <- expand.grid(X = (1:n) + (img-1)*n, Y = 1:n)
}
data_combined <- do.call(rbind, grid_list)
colnames(data_combined) <- c("X_Coordinate", "Y_Coordinate")
N_points <- nrow(data_combined)

# Distance matrix
m <- computeDistanceMatrix(as.matrix(data_combined))

# Prior kernel visualization
p.etasq <- rexp(N_points, 2)[1:20]
p.rhosq <- rexp(N_points, 0.5)[1:20]
pdf(file.path(out_folder, "kernel_prior.pdf"))
plot(NULL, xlim = c(0, max(m)/3), ylim = c(0, 1),
     xlab = "pixel distance", ylab = "covariance", main = "Prior simulation")
for(i in 1:20)
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2), add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
dev.off()

# GP simulation
beta <- 5
etasq <- 2
rho <- sqrt(0.5)
K <- etasq * exp(-0.5 * ((m / rho)^2)) + diag(1e-9, N_points)
sim_gp <- MASS::mvrnorm(1, mu = rep(0, N_points), Sigma = K)
sim_y <- rnorm(N_points, mean = sim_gp + beta * x_all, sd = 1)

# MALDI images
pdf(file.path(out_folder, "simulated_maldi_combined.pdf"))
par(mfrow = c(ceiling(sqrt(n_images)), ceiling(sqrt(n_images))))
for(img in 1:n_images) {
  start_idx <- (img-1) * (n * n) + 1
  end_idx <- img * (n * n)
  image(matrix(sim_y[start_idx:end_idx], nrow = n, ncol = n, byrow = TRUE),
        col = viridis, main = sprintf("Simulated MALDI - Image %d", img),
        xlab = "", ylab = "", axes = FALSE)
}
dev.off()

# Data preparation
dat_list <- list(y = sim_y, x = x_all, S = 1:N_points, Dmat = m, N = N_points)

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
  data = dat_list, chains = 4, cores = 4, iter = 1000
)

# Print summary
precis(GP)
post <- extract.samples(GP)

# Kernel plot
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
