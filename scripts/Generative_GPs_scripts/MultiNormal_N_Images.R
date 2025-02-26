set.seed(20241028)
library(rethinking)
library(plot.matrix)
library(viridis)
library(MASS)

# Function to simulate a vessel mask matrix
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

# -------------------------
# Set N of Images
# -------------------------
N_img <- 5        # Runs with 5 --> 1200 sec needed for 250 iterations
n <- 20           # Matrix dimensions (n x n)

circles_list <- vector("list", N_img)
for(i in 1:N_img) {
  circles_list[[i]] <- list(
    list(center = c(sample(5:15, 1), sample(5:15, 1)), radius = runif(1, 1, 3)),
    list(center = c(sample(5:15, 1), sample(5:15, 1)), radius = runif(1, 1, 3))
  )
}

# Generate vessel mask matrices for each image
mats <- vector("list", N_img)
for(i in 1:N_img) {
  mats[[i]] <- test_function(n, circles_list[[i]])
}

# Transform matrices into vectors
xs <- vector("list", N_img)
for(i in 1:N_img) {
  xs[[i]] <- as.vector(t(mats[[i]]))
}

# Create coordinate grids for each image
grids <- vector("list", N_img)
for(i in 1:N_img) {
  grids[[i]] <- expand.grid(X = (1:n) + (i-1)*n, Y = 1:n)
}

# Combine all grids to form the global dataset
data_combined <- do.call(rbind, grids)
colnames(data_combined) <- c("X_Coordinate", "Y_Coordinate")
N_points <- nrow(data_combined)

# Compute the global distance matrix
m_global <- as.matrix(dist(data_combined, method = "euclidean"))

# Compute individual distance matrices for each image
Dmats <- vector("list", N_img)
for(i in 1:N_img) {
  Dmats[[i]] <- as.matrix(dist(grids[[i]], method = "euclidean"))
}

# -------------------------
# SIM
# -------------------------
beta <- 5
etasq <- 2
rho <- sqrt(0.5)

# Generate covariance matrices for each image
Ks <- vector("list", N_img)
for(i in 1:N_img) {
  Ks[[i]] <- etasq * exp(-0.5 * ((Dmats[[i]] / rho)^2)) + diag(1e-9, n*n)
}

# Sample from the GP prior
sim_gp <- vector("list", N_img)
for(i in 1:N_img) {
  sim_gp[[i]] <- MASS::mvrnorm(1, mu = rep(0, n*n), Sigma = Ks[[i]])
}

# Generate observed values
sim_y <- vector("list", N_img)
for(i in 1:N_img) {
  sim_y[[i]] <- rnorm(n*n, mean = sim_gp[[i]] + beta * xs[[i]], sd = 1)
}

# -------------------------
# Dati prep
# -------------------------
dat_list <- list(N = n*n)  # number of points per image
for(i in 1:N_img) {
  dat_list[[ paste0("y", i) ]] <- sim_y[[i]]
  dat_list[[ paste0("x", i) ]] <- xs[[i]]
  dat_list[[ paste0("Dmat", i) ]] <- Dmats[[i]]
}

# -------------------------
#Covariance Matrices
# -------------------------
model_code <- "alist(\n"
for(i in 1:N_img) {
  model_code <- paste0(model_code,
                       "  y", i, " ~ multi_normal(mu", i, ", K", i, "),\n",
                       "  mu", i, " <- a + b * x", i, ",\n",
                       "  matrix[N, N]:K", i, " <- cov_GPL2(Dmat", i, ", etasq, rho, 0.01),\n")
}
model_code <- paste0(model_code,
                     "  a ~ normal(0, 1),\n",
                     "  b ~ normal(0, 0.5),\n",
                     "  etasq ~ dexp(2),\n",
                     "  rho ~ dexp(0.5)\n",
                     ")")

cat(model_code)

model_list <- eval(parse(text = model_code))

# -------------------------
#  FIT THE MODEL
# -------------------------
GP_N <- ulam(model_list, data = dat_list, chains = 4, cores = 4, iter = 250)

print(precis(GP_N))
post <- extract.samples(GP_N)

# -------------------------
# PLOT KERNEL COMPARISON
# -------------------------
# Plot prior curves, the actual kernel, and estimated kernel curves 
par(mfrow = c(1, 1))
plot(NULL, xlim = c(0, max(m_global)/3), ylim = c(0, 2),
     xlab = "pixel distance", ylab = "covariance",
     main = "Prior, Actual, and Estimated Kernel")
# Plot priors
for(i in 1:20)
  curve(p.etasq[i] * exp(-p.rhosq[i] * x^2),
        add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
# Plot actual kernel
curve(etasq * exp(-0.5 * rho * x^2), add = TRUE, lwd = 4)
# Plot estimated kernels (using 20 posterior samples)
for(i in 1:20) {
  curve(post$etasq[i] * exp(-post$rho[i] * x),
        add = TRUE, col = col.alpha(4, 0.3), lwd = 6)
}
legend("topright", lwd = 4, col = c(2,1,4),
       legend = c("Prior", "Actual", "Estimated"))
