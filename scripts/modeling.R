library(rethinking)

# Read data
db <- read.csv('./input/db.csv', stringsAsFactors = FALSE)

# Process data
db$model <- as.factor(db$model)
db$treatment <- as.factor(db$treatment)
db$min_log <- log(db$axis_minor_length)


# Prepare data for modeling
d <- list(
  axis_length = standardize(db$min_log),
  mod = as.numeric(db$model), # 1: "LM02", 2: "LM04"
  rx = as.numeric(db$treatment) # 1: "CTR", 2: "DTX", 3: "ED1", 4: "ED2"
)

#########################  1  ##############################

# Model 1: Exponential distribution with full interaction model
# This allows for different effects of ERI for each model
m <-ulam(
  alist(
    axis_length ~ dnorm(mu, sigma),
    mu <- a[mod] + b[rx],
    a[mod] ~ dnorm(a_bar, sigma_a),
    b[rx] ~ dnorm(b_bar, sigma_b),
    c(a_bar, b_bar) ~ dnorm(0,0.7),
    c(sigma,sigma_a,sigma_b) ~ dexp(1)
  ), data=d, chains=4, cores=4, iter= 1000)


# 
# mnc <- ulam(
#   alist(
#     axis_length ~ dnorm(mu, zeta),
#     mu <- a[mod] + b[rx],
#     # define effects using other parameters
#     save> vector[50633]:a <<- abar + za*sigma,
#     save> vector[50633]:b <<- bbar + zb*tau,
#     # z-scored effects
#     vector[50633]:za ~ normal(0,1),
#     vector[50633]:zb ~ normal(0,1),
#     # ye olde hyper-priors
#     c(abar,bbar) ~ normal(0,1),
#     c(sigma,tau, zeta) ~ exponential(1)
#   ) , data=d , chains=4 , cores=4 , iter = 500)

# Print model results
dashboard(m)
precis(m)

# Plot posterior distributions
post <- extract.samples(m)


p_link <- function( mod=1 , rx=1 ) {
  x <- with( post , a[,mod] + b[,rx] )
  return( x )
}


dest <- function (x){
   x <- x * 0.378
   x <- x + 4.44
   return(x)
}

post_mu <- p_link(d$mod,d$rx)
sim_min <- rnorm(length(d$axis_length), post_mu[1,], 1 )
sim_min <- dest(sim_min)
plot(NULL, xlim = c(0,200), ylim = c(0,0.03), 
     main = 'Posterior Predictive Check', 
     xlab = 'Minor Axis Size',
     ylab = 'Density')
dens(db$axis_minor_length, lwd = 3, add = TRUE)
dens(exp(sim_min), lwd = 3, add = TRUE, col = 2)


# no model
p_link <- function(  rx=1 ) {
  x <- with( post , rnorm(length(post$a_bar), b[,rx], sigma_b) )
  x <- dest(x)
  return( x )
}

p_raw_bsim <- sapply(1:4, function(i) p_link(i))
#conttrast
p_raw_bsim[,1]


png('output/contrastRXnoModel.png', 
    width = 3000, height = 2000, res = 300)
par(mai = c(1.1, 0.8, 0.8, 0.4) + 0.02)
col <- c("#332288","#AA4499","#117733")
s_factor <- 20 # scaling factor for graphics
new.lab = c("DTX", "ERI1", "ERI2")
plot(NULL, ylim=c(-0.8,1), xlim = c(0.5,3.5),
     xlab = '', ylab = 'Probability Density', xaxt = 'null', 
     main = 'Treatment Contrast vs CTRL')
abline(h = 0)
for(i in 1:3) {
  obj <- p_raw_bsim[,i+1] - p_raw_bsim[,1]
  y <- density(obj)$x
  x <- density(obj)$y
  polygon(i + x/s_factor, y, col = scales::alpha(col[[i]],0.6), border = FALSE)
  lines(i + x/s_factor, y, lwd = 1)
  polygon(i - x/s_factor, y, col = scales::alpha(col[[i]],0.6), lwd = 2, border = FALSE)
  lines(i - x/s_factor, y, lwd = 1)
}
axis(1, at = 1:3, labels = FALSE)
text(x = 1:3+0.05, -0.92,
     labels = new.lab,  
     srt=40,  adj=1,    xpd=TRUE)
dev.off()
