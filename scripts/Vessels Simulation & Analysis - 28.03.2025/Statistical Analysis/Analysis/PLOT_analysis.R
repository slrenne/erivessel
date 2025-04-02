library(rethinking)

# Read data - using a different name than 'data' to avoid conflicts
my_data <- read.csv("input/db.csv", stringsAsFactors = FALSE)


i <- my_data$treatment %in% c( "ED2", "ED1" )
my_data$model <- as.factor(my_data$model)
my_data$RxE <- i + 1L # 2 ERI 1 else
#i <- my_data$model == "LM02"

png("output/minAxSize_density.png")
plot(NULL, xlim = c(0,300), ylim= c(0,0.03), xaxt = "n", yaxt = "n", xlab = "Size", ylab = "Density", main = "Minor Axis")
for(i in 1:2){
  for(j in 1:2){
    idx_i <- my_data$model == levels(my_data$model)[i]
    idx_j <- my_data$RxE == j
    idx <- idx_i & idx_j
    d <- density(my_data$axis_minor_length[idx], adjust = 0.5)
    lines(d, add = TRUE, lwd = 2, col = i, lty = j)
  }
}
legend("topright", 
       lty = c(1,2,1,2), 
       lwd = rep(2, times = 4), 
       col = c(1,1,2,2), 
       legend = c("LM02 else", "LM02 ERI", "LM04 else", "LM04 ERi"))
dev.off()

###############################
# modeling##

dat <- list(
  a = my_data$axis_minor_length, # target minor axis
     )
m_e <- ulam(
  alist(
    a ~ dexp(m),
    m <- 
    m ~ dexp(1)
    ), data=dat , chains=4 , cores=4)

precis(m_e)
