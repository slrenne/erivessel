set.seed(20241028)
library(rethinking)
library(plot.matrix)
library(viridis)
library(MASS)
library(parallel)  #  parallel library

# Set up parallel cluster for Windows
num_cores <- detectCores() - 1  # Use all but one core 
cl <- makeCluster(num_cores)    # Create a parallel cluster

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
N_img <- 10       # Runs with 10 --> 30 min for 600 iterations
n <- 20           # Matrix dimensions (n x n)

# Generate random circles for each image
circles_list <- vector("list", N_img)
for(i in 1:N_img) {
  circles_list[[i]] <- list(
    list(center = c(sample(5:15, 1), sample(5:15, 1)), radius = runif(1, 1, 3)),
    list(center = c(sample(5:15, 1), sample(5:15, 1)), radius = runif(1, 1, 3))
  )
}

# Export necessary variables to the cluster
clusterExport(cl, c("test_function", "n", "circles_list"))

# Generate vessel mask matrices for each image in parallel
mats <- parLapply(cl, 1:N_img, function(i) test_function(n, circles_list[[i]]))

# Export mats to the cluster
clusterExport(cl, "mats")

# Transform matrices into vectors in parallel
xs <- parLapply(cl, 1:N_img, function(i) as.vector(t(mats[[i]])))

# Create coordinate grids for each image in parallel
grids <- parLapply(cl, 1:N_img, function(i) expand.grid(X = 1:n, Y = 1:n))

# Export grids to the cluster
clusterExport(cl, "grids")

# Compute individual distance matrices for each image in parallel
Dmats <- parLapply(cl, 1:N_img, function(i) as.matrix(dist(grids[[i]], method = "euclidean")))

# -------------------------
# SIM
# -------------------------
beta <- 5
etasq <- 2
rho <- sqrt(0.5)

# Export simulation parameters to the cluster
clusterExport(cl, c("beta", "etasq", "rho", "Dmats"))

# Generate covariance matrices for each image in parallel
Ks <- parLapply(cl, 1:N_img, function(i) {
  etasq * exp(-0.5 * ((Dmats[[i]] / rho)^2)) + diag(1e-9, n*n)
})

# Export Ks to the cluster
clusterExport(cl, "Ks")

# Sample from the GP prior for each image independently in parallel
sim_gp <- parLapply(cl, 1:N_img, function(i) {
  MASS::mvrnorm(1, mu = rep(0, n*n), Sigma = Ks[[i]])
})

# Export sim_gp to the cluster
clusterExport(cl, c("xs", "sim_gp"))

# Generate observed values in parallel
sim_y <- parLapply(cl, 1:N_img, function(i) {
  rnorm(n*n, mean = sim_gp[[i]] + beta * xs[[i]], sd = 1)
})

# -------------------------
# Data prep
# -------------------------
dat_list <- list(N = n*n)  # number of points per image
for(i in 1:N_img) {
  dat_list[[ paste0("y", i) ]] <- sim_y[[i]]
  dat_list[[ paste0("x", i) ]] <- xs[[i]]
  dat_list[[ paste0("Dmat", i) ]] <- Dmats[[i]]
}

# -------------------------
# Covariance Matrices
# -------------------------

model_code <- "alist(\n"
for(i in 1:N_img) {
  model_code <- paste0(model_code,
                       "  y", i, " ~ multi_normal(mu", i, ", K", i, "),\n",
                       "  mu", i, " <- a + b * x", i, ",\n",
                       "  matrix[N, N]:K", i, " <- etasq * exp(-0.5 * square(Dmat", i, " / rho)) + diag_matrix(rep_vector(0.01, N)),\n")
}
model_code <- paste0(model_code,
                     "  a ~ normal(0, 1),\n",
                     "  b ~ normal(0, 0.5),\n",
                     "  etasq ~ exponential(2),\n",
                     "  rho ~ exponential(0.5)\n",
                     ")")

cat(model_code)

model_list <- eval(parse(text = model_code))

# -------------------------
#  FIT THE MODEL
# -------------------------
# Using the parallel cores for model fitting
GP_N <- ulam(model_list, data = dat_list, chains = 4, cores = num_cores, iter = 600, warmup = 150)

print(precis(GP_N))
post <- extract.samples(GP_N)

# -------------------------
# PLOT KERNEL COMPARISON
# -------------------------
# sample from the priors
set.seed(08062002) #re set seed otherwise prior curves keep changing 

p.etasq <- rexp(n, rate = 0.5)
p.rhosq <- rexp(n, rate = 0.5)

# Plot prior curves, the actual kernel, and estimated kernel curves 
par(mfrow = c(1, 1))
plot(NULL, xlim = c(0, max(Dmats[[1]])/3), ylim = c(0, 10),
     xlab = "pixel distance", ylab = "covariance",
     main = "Prior, Actual, and Estimated Kernel")
# Plot priors
for(i in 1:20)
  curve(p.etasq[i] * exp(-0.5 * (x/p.rhosq[i])^2),
        add = TRUE, lwd = 6, col = col.alpha(2, 0.5))
# Plot actual kernel
curve(etasq * exp(-0.5 * (x/rho)^2), add = TRUE, lwd = 4)
# Plot estimated kernels (using 20 posterior samples)
for(i in 1:20) {
  curve(post$etasq[i] * exp(-0.5 * (x/post$rho[i])^2),
        add = TRUE, col = col.alpha(4, 0.3), lwd = 6)
}
legend("topright", lwd = 4, col = c(2,1,4),
       legend = c("Prior", "Actual", "Estimated"))

# -------------------------
# STOP THE CLUSTER (IMPORTANT)
# -------------------------
stopCluster(cl)
