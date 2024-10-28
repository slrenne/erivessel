set.seed(20241028)

# Simulate spatially autocorrelated Poisson data
library(rethinking)
library(plot.matrix)
library(viridis)


# Set matrix dimensions
n <- 20
mat <- matrix(0, nrow = n, ncol = n)
X_coordinate <- 1:n
Y_coordinate <- 1:n
data <- expand.grid(X = X_coordinate, Y = Y_coordinate)
colnames(data) <- c('X_Coordinate', 'Y_Coordinate')
N_points <- nrow(data)

# Define center and radius for two circles
center1 <- c(7, 7)  # Center of the first circle
radius1 <- 1        # Radius of the first circle
center2 <- c(14, 14)  # Center of the second circle
radius2 <- 2         # Radius of the second circle

# Fill the matrix to create circles
for (i in 1:n) {
  for (j in 1:n) {
    # Calculate distance from the first circle's center
    distance1 <- sqrt((i - center1[1])^2 + (j - center1[2])^2)
    # Calculate distance from the second circle's center
    distance2 <- sqrt((i - center2[1])^2 + (j - center2[2])^2)
    
    # Set matrix values based on whether the point is within each circle
    if (distance1 <= radius1) mat[i, j] <- 1
    if (distance2 <= radius2) mat[i, j] <- 1
  }
}

# Plot the matrix with colors
pdf("./output/simulated_vessel.pdf")
plot(mat, col = viridis, key = NULL, main = "Simulated vessels",
     xlab="", ylab="", axis.col=NULL, axis.row=NULL)
dev.off()

# transforming into the vector
x <- vector()
for (i in 1:n) {
  for (j in 1:n) {
    x[ (i - 1) * 20 + j] <- mat[i,j]
  }}
plot(x)

# calculate the distance matrix
m <- as.matrix(dist(data, method = "euclidean"))

# sim priors for distance model
p.etasq <- rexp(n,2)
p.rhosq <- rexp(n,0.5)

pdf("./output/kernel_prior.pdf")
plot( NULL , xlim=c(0,max(m)/3) , ylim=c(0,1) , xlab="pixel distance" , ylab="covariance" , main = "Prior simulation")
for ( i in 1:n )curve( p.etasq[i]*exp(-p.rhosq[i]*x^2) , add=TRUE , lwd=6 , col=col.alpha(2,0.5) )
dev.off()

# setting the parameters
beta <- 5
etasq <- 2
rho <- sqrt(0.5)


# Generate Gaussian Process
K <- etasq  * exp(-0.5 * ((m / rho) ^ 2)) + diag(1e-9, N_points) # @Omer why 0.5?


# Sample from MVNORM
sim_gp <- MASS::mvrnorm(1, mu = rep(0, N_points), Sigma = K)

sim_y <- rnorm(N_points, sim_gp + beta * x)
pdf("./output/simulated_maldi.pdf")
plot(matrix(sim_y, ncol = n, byrow = TRUE), col = viridis, 
     key = NULL, main = "Simulated MALDI",
     xlab="", ylab="", axis.col=NULL, axis.row=NULL)
dev.off()

# gaussian process
dat_list <- list(
  y = sim_y,
  x = x,
  S = 1:N_points,
  Dmat = m , 
  N = N_points)


GP <- ulam(
  alist(
    y ~ multi_normal( mu , K), 
    mu <-  a + b* x,
    matrix[N,N]:K <- cov_GPL2(Dmat,etasq,rho,0.01),
    a ~ normal(0,1),
    b ~ normal(0,0.5),
    etasq ~ dexp(2),
    rho ~ dexp(0.5)
  ), data=dat_list , chains=4 , cores=4 , iter=1000 )

precis(GP)

post <- extract.samples(GP)

# plotting the kernel 
pdf("./output/kernel_prior_actual_post.pdf")
plot( NULL , xlim=c(0,max(m)/3) , ylim=c(0,2), 
      xlab="pixel distance" , ylab="covariance" , 
      main = "Prior, Actual, and Estimated Kernel")


for ( i in 1:20 )curve( p.etasq[i]*exp(-p.rhosq[i]*x^2) , add=TRUE , lwd=6 , col=col.alpha(2,0.5) )

x <- NULL
curve( etasq * exp(-0.5 * rho * x^2) , add=TRUE , lwd=4 )

for ( i in 1:20 ) {  
  curve( post$etasq[i]*exp(-post$rho[i]*x) , 
         add=TRUE , col=col.alpha(4,0.3) , lwd= 6 )
} 
legend("topright", lwd = 4, col = c(2,1,4), legend = c("Prior","Actual","Estimated"))
dev.off()