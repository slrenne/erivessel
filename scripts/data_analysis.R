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

# count the vessels
v_count <- data.frame(table(db$Scan))
colnames(v_count) <- c("Scan","count")
key <- read.delim('input/scan_model_mouse_rx.csv', sep = ";")
cdb <- merge(v_count, key, by = "Scan")
cdb$eribulin <- ifelse(cdb$treatment %in% c("ED1","ED2"), 
                      1L, # had eribulin
                      0L) # did not 

for (i in 1:5) {cdb[,i+2] <- as.factor(cdb[,i+2])}

dat <- list(
  V = standardize(cdb$count),
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

pdf("./output/totvasc.pdf")
plot( NULL , 
      xlab="Eribulin" , 
      xlim=c(0.5,2.5), 
      ylim=c(-1.1,1) , 
      ylab="vessel count",
      xaxt = "n",
      yaxt = "n",
      main = "Absolute number of vessels")
ytick <- seq(-1.1, 1, length = 5)
axis(2, at = ytick, round(de_standard(ytick, dat$V)))
axis(1, at = 1:2, c("no","yes"))
for ( i in 1:2 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(2,0.5) )
points( 1:2 , apply(p,2,mean), lwd=3 , col=2 , pch = 16)
dev.off()

pdf("./output/totvasc_contrast.pdf")
mu_cont <- de_standard(p[,2], dat$V) - de_standard(p[,1], dat$V)
dens(mu_cont, lwd = 3, 
     xlab = "number of vessels gained with eribulin",
     main  = "Posterior mean contrast")
dev.off()

