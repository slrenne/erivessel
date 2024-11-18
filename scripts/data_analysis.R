# set seed
set.seed(240930)

library(rethinking)
de_standard <- function(x, list_norm){ 
  x <- x * attr(list_norm,"scaled:scale")
  x <- x + attr(list_norm,"scaled:center")
  x
  }

# load the db
db <- read.csv("input/db.csv")
# add a colun for Eribuin rx
db$eribulin <- ifelse(db$treatment %in% c("ED1","ED2"), 
                      1L, # had eribulin
                      0L) # did not 

db$eribulin2 <- ifelse(db$treatment == "ED1",
                      1L,ifelse(db$treatment == "ED2",
                                2L, # had eribulin
                                0L)) # did not 


# count the vessels
v_count <- data.frame(table(db$Scan))
colnames(v_count) <- c("Scan","count")
key <- read.delim('input/scan_model_mouse_rx_area.csv', sep = ";")
key$area.px.2. <- as.integer(key$area.px.2.)
cdb <- merge(v_count, key, by = "Scan")
cdb$eribulin <- ifelse(cdb$treatment %in% c("ED1","ED2"), 
                      1L, # had eribulin
                      0L) # did not 

for (i in c(1:4,6)) {cdb[,i+2] <- as.factor(cdb[,i+2])}

# modeling the total number of vessels 

dat <- list(
  V = standardize(cdb$count / cdb$area.px.2.),
  E = as.integer(cdb$eribulin) )
# centered version of the model
# m <- ulam(
#   alist(
#     V ~ dnorm(mu,sigma),
#     mu <- a[E] + b[E],
#     vector[2]:a ~ normal(abar,sigma),
#     vector[2]:b ~ normal(bbar,tau),
#     c(abar,bbar) ~ normal(0,1),
#     c(sigma,tau) ~ exponential(1)
#   ), data=dat  , chains=4 , cores=4)
# non centered
m <- ulam(
  alist(
    V ~ dnorm(mu,sigma),
    mu <- a[E] + b[E],
    # define effects using other parameters
    save> vector[2]:a <<- abar + za*sigma,
    save> vector[2]:b <<- bbar + zb*tau,
    # z-scored effects
    vector[2]:za ~ normal(0,1),
    vector[2]:zb ~ normal(0,1),
    # ye olde hyper-priors
    c(abar,bbar) ~ normal(0,1),
    c(sigma,tau) ~ exponential(1)
  ), data=dat  , chains=4 , cores=4)

p <- link(m, data = list(E = 1:2))

mu_cont <- p[,2] - p[,1]



pdf("./output/totvasc.pdf")
# Set up a 1x2 layout, where the first plot is on the left and the second is on the right
l.m <- matrix(c(1, 2), nrow = 1, ncol = 2) # 1 row, 2 columns
layout(mat = l.m, widths = c(10, 4))          # adjust widths as desired

# First plot (left side)
plot(NULL, 
     xlab = "Eribulin", 
     xlim = c(0.5, 2.5), 
     ylim = c(-1.1, 1), 
     ylab = "Normalized Vessels' Density",
     xaxt = "n",
     main = "Vessels' Density")
ytick <- seq(from = -1.5, to = 1.5, length = 5)
axis(1, at = 1:2, labels = c("no", "yes"))
for (i in 1:2) lines(c(i, i), PI(p[, i]), lwd = 8, col = col.alpha(2, 0.5))
points(1:2, apply(p, 2, mean), lwd = 3, col = 2, pch = 16)

# Second plot (right side), rotated by 90 degrees
par(mar = c(5, 5, 4, 2) + 0.1) # adjust margins if needed
# Rotate the second plot by setting `par` to flip coordinates
par(srt = 90) # angle text
a <- density(mu_cont,adjust = 0.5)

plot(a$y, a$x, type = 'l', lwd = 4,
     xlab = "Density", 
     ylab = "Expected Increase in Vessels'\n Density with Eribulin",
     main  = "",
     xlim = c(0, 0.9),
     ylim = c(-2.2, 2))
abline(h = 0, lty = 2)


# Close any open graphics device
dev.off()





# modeling minor axis 

dat <- list(
  mL = standardize(db$axis_minor_length),
  E = db$eribulin + 1L )

m <- quap(
  alist(
    mL ~ dnorm(mu,sigma),
    mu <- a[E] + b[E],
    a[E] ~ dnorm(0,1),
    b[E] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=dat)
# # centered version of the model
# m <- ulam(
#   alist(
#     mL ~ dnorm(mu,sigma),
#     mu <- a[E] + b[E],
#     vector[2]:a ~ normal(abar,sigma),
#     vector[2]:b ~ normal(bbar,tau),
#     c(abar,bbar) ~ normal(0,1),
#     c(sigma,tau) ~ exponential(1)
#   ), data=dat  , chains=4 , cores=4)
# # non centered
# m <- ulam(
#   alist(
#     mL ~ dnorm(mu,sigma),
#     mu <- a[E] + b[E],
#     # define effects using other parameters
#     save> vector[2]:a <<- abar + za*sigma,
#     save> vector[2]:b <<- bbar + zb*tau,
#     # z-scored effects
#     vector[2]:za ~ normal(0,1),
#     vector[2]:zb ~ normal(0,1),
#     # ye olde hyper-priors
#     c(abar,bbar) ~ normal(0,1),
#     c(sigma,tau) ~ exponential(1)
#   ), data=dat  , chains=4 , cores=4)

p <- link(m, data = list(E = 1:2))

pdf("./output/miL.pdf")
plot( NULL , 
      xlab="Eribulin" , 
      xlim=c(0.5,2.5), 
      ylim=c(min(p),max(p)) , 
      ylab= expression("length (in px)"),
      xaxt = "n",
      yaxt = "n",
      main = "Minor length axis")
ytick <- seq(min(p),max(p), length = 5)
axis(2, at = ytick, round(de_standard(ytick, dat$mL)))
axis(1, at = 1:2, c("no","yes"))
for ( i in 1:2 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(2,0.5) )
points( 1:2 , apply(p,2,mean), lwd=3 , col=2 , pch = 16)
dev.off()

pdf("./output/miL_contrast.pdf")
mu_cont <- de_standard(p[,2], dat$mL) - de_standard(p[,1], dat$mL)
dens(mu_cont, lwd = 3, 
     xlab =  expression("length gained with eribulin (in px)"),
     main  = "Posterior mean contrast minor axis", show.zero = TRUE)
dev.off()


# modeling Major axis 

dat <- list(
  Ec = standardize(db$axis_major_length),
  E = db$eribulin + 1L )

m <- quap(
  alist(
    Ec ~ dnorm(mu,sigma),
    mu <- a[E] + b[E],
    a[E] ~ dnorm(0,1),
    b[E] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=dat)

p <- link(m, data = list(E = 1:2))

pdf("./output/MaL.pdf")
plot( NULL , 
      xlab="Eribulin" , 
      xlim=c(0.5,2.5), 
      ylim=c(min(p),max(p)) , 
      ylab= expression("length  (in px)"),
      xaxt = "n",
      yaxt = "n",
      main = "Major length axis")
ytick <- seq(min(p),max(p), length = 5)
axis(2, at = ytick, round(de_standard(ytick, dat$Ec)))
axis(1, at = 1:2, c("no","yes"))
for ( i in 1:2 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(2,0.5) )
points( 1:2 , apply(p,2,mean), lwd=3 , col=2 , pch = 16)
dev.off()

pdf("./output/MaL_contrast.pdf")
mu_cont <- de_standard(p[,2], dat$Ec) - de_standard(p[,1], dat$Ec)
dens(mu_cont, lwd = 3, 
     xlab =  expression("length gained with eribulin (in px)"),
     main  = "Posterior mean contrast major axis", show.zero = TRUE)
dev.off()

# modeling eccentricity

dat <- list(
  Ec = standardize(db$eccentricity),
  E = db$eribulin + 1L )

m <- quap(
  alist(
    Ec ~ dnorm(mu,sigma),
    mu <- a[E] + b[E],
    a[E] ~ dnorm(0,1),
    b[E] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=dat)

p <- link(m, data = list(E = 1:2))

pdf("./output/Ec.pdf")
plot( NULL , 
      xlab="Eribulin" , 
      xlim=c(0.5,2.5), 
      ylim=c(min(p),max(p)) , 
      ylab= "Eccentricity",
      xaxt = "n",
      yaxt = "n",
      main = "Eccentricity")
ytick <- seq(min(p),max(p), length = 5)
axis(2, at = ytick, round(de_standard(ytick, dat$Ec)))
axis(1, at = 1:2, c("no","yes"))
for ( i in 1:2 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(2,0.5) )
points( 1:2 , apply(p,2,mean), lwd=3 , col=2 , pch = 16)
dev.off()

pdf("./output/Ec_contrast.pdf")
mu_cont <- de_standard(p[,2], dat$Ec) - de_standard(p[,1], dat$Ec)
dens(mu_cont, lwd = 3, 
     xlab =  "Eccentricity gained with eribulin",
     main  = "Posterior mean contrast minor axis", show.zero = TRUE)
dev.off()


# modelingvessels area 

dat <- list(
  A = standardize(db$area),
  E = db$eribulin + 1L )

m <- quap(
  alist(
    A ~ dnorm(mu,sigma),
    mu <- a[E] + b[E],
    a[E] ~ dnorm(0,1),
    b[E] ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=dat)

p <- link(m, data = list(E = 1:2))

pdf("./output/A.pdf")
plot( NULL , 
      xlab="Eribulin" , 
      xlim=c(0.5,2.5), 
      ylim=c(min(p),max(p)) , 
      ylab= expression("Area  (in px" ^2*")"),
      xaxt = "n",
      yaxt = "n",
      main = "Vessels' Area")
ytick <- seq(min(p),max(p), length = 5)
axis(2, at = ytick, round(de_standard(ytick, dat$A)))
axis(1, at = 1:2, c("no","yes"))
for ( i in 1:2 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(2,0.5) )
points( 1:2 , apply(p,2,mean), lwd=3 , col=2 , pch = 16)
dev.off()

pdf("./output/A_contrast.pdf")
mu_cont <- de_standard(p[,2], dat$A) - de_standard(p[,1], dat$A)
dens(mu_cont, lwd = 3, 
     xlab =  expression("Area gained with eribulin (in px" ^2*")"),
     main  = "Posterior mean contrast area", show.zero = TRUE)
dev.off()

