
while(dev.cur() > 1) dev.off()

set.seed(20241028)
library(MASS)
library(viridis)
library(parallel)
library(ggplot2)
library(tibble)
library(bitops) 
library(gghilbertstrings)  
library(pROC) 
library(plot.matrix)
library(viridis)

# 1. Simulation of Vessel and MALDI Images

test_function <- function(n, circles) {
  mat <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      for(circle in circles) {
        center <- circle$center
        radius <- circle$radius
        if (sqrt((i - center[1])^2 + (j - center[2])^2) <= radius)
          mat[i, j] <- 1
      }
    }
  }
  return(mat)
}

N_img_per_group <- 4    
n <- 20       

circles_list <- vector("list", 3 * N_img_per_group) 
vessel_images <- vector("list", 3 * N_img_per_group)

for(i in 1:(3 * N_img_per_group)) {
  circles_list[[i]] <- list(
    list(center = c(sample(5:15, 1), sample(5:15, 1)), radius = runif(1, 1, 3)),
    list(center = c(sample(5:15, 1), sample(5:15, 1)), radius = runif(1, 1, 3))
  )
  vessel_images[[i]] <- test_function(n, circles_list[[i]])
}

# Parameters for the different groups
beta_strong <- 5     # Strong correlation (original)
beta_weak <- 1       # Weak correlation
beta_none <- 0       # No correlation (pure noise)
etasq <- 2
rho <- sqrt(0.5)

grids <- expand.grid(X = 1:n, Y = 1:n)
Dmat <- as.matrix(dist(grids, method = "euclidean"))
K <- etasq * exp(-0.5 * ((Dmat / rho)^2)) + diag(1e-9, n*n)

maldi_images <- vector("list", 3 * N_img_per_group)
image_groups <- rep(1:3, each = N_img_per_group) # Group labels

for(i in 1:N_img_per_group) {
  nn <- n*n  # Fix: Define nn
  gp_i <- MASS::mvrnorm(1, mu = rep(0, nn), Sigma = K)
  maldi_vec <- rnorm(nn, mean = gp_i + beta_strong * as.vector(t(vessel_images[[i]])), sd = 1)
  maldi_images[[i]] <- matrix(maldi_vec, nrow = n, ncol = n, byrow = TRUE)
}

# Group 2: Weak correlation
for(i in 1:N_img_per_group) {
  nn <- n*n
  gp_i <- MASS::mvrnorm(1, mu = rep(0, nn), Sigma = K)
  maldi_vec <- rnorm(nn, mean = gp_i + beta_weak * as.vector(t(vessel_images[[i + N_img_per_group]])), sd = 1)
  #add a diffusion param to posterior sim [to add diffusion to Hilbertreg]
  maldi_images[[i + N_img_per_group]] <- matrix(maldi_vec, nrow = n, ncol = n, byrow = TRUE)
}

# Group 3: No correlation (random noise)
for(i in 1:N_img_per_group) {
  nn <- n*n
  gp_i <- MASS::mvrnorm(1, mu = rep(0, nn), Sigma = K)
  maldi_vec <- rnorm(nn, mean = gp_i + beta_none * as.vector(t(vessel_images[[i + 2*N_img_per_group]])), sd = 1)
  maldi_images[[i + 2*N_img_per_group]] <- matrix(maldi_vec, nrow = n, ncol = n, byrow = TRUE)
}


###################################
plot(maldi_images[[1]])
pdf("output/simulated_pairs_strongWeakNoise.pdf", height = 12, width = 8)
par(mfrow=c(6,4))
for(j in c(0,4,8)){
  for(i in j+1:4) { plot(vessel_images[[i]], 
                     col = magma, key = NULL, main = "Simulated vessels",
     xlab="", ylab="", axis.col=NULL, axis.row=NULL)}
  for(i in j+1:4) {plot(maldi_images[[i]], 
                    col = viridis, key = NULL, main = "Simulated maldi",
     xlab="", ylab="", axis.col=NULL, axis.row=NULL)}
}
dev.off()

#############################################
plot_hilbert_image <- function(image_matrix, title_text) {
  vec <- as.vector(t(image_matrix))
  
  df <- tibble(val = 1:(n*n), value = vec)
  
  p <- gghilbertplot(df, val , size = 5, 
                     color = value, add_curve = TRUE) +
    scale_color_viridis() +
    ggtitle(title_text) +
    theme_minimal()
  return(p)
}

display_image_pair <- function(index, save_plots = FALSE) {
  group_name <- c("Strong Correlation", "Weak Correlation", "No Correlation")[image_groups[index]]
  
  vessel_plot <- plot_hilbert_image(vessel_images[[index]], 
                                    paste("Vessel Image -", group_name, "- Image", index))
  
  maldi_plot <- plot_hilbert_image(maldi_images[[index]], 
                                   paste("MALDI Image -", group_name, "- Image", index))
  
  if(!save_plots) {
    print(vessel_plot)
    print(maldi_plot)
  } else {
    # Save plots if required
    filename_base <- paste0("Image_", index, "_", gsub(" ", "_", group_name))
    ggsave(paste0(filename_base, "_vessel.png"), vessel_plot)
    ggsave(paste0(filename_base, "_maldi.png"), maldi_plot)
  }
}

for(i in 1:(3 * N_img_per_group)) {
  display_image_pair(i)
}



xy2d <- function(n_grid, x, y) {
  d <- 0
  s <- n_grid / 2
  while(s >= 1) {
    rx <- ifelse(bitwAnd(x, s) != 0, 1, 0)
    ry <- ifelse(bitwAnd(y, s) != 0, 1, 0)
    d <- d + s * s * ((3 * rx) ^ ry)
    if (ry == 0) {
      if(rx == 1) {
        x <- n_grid - 1 - x
        y <- n_grid - 1 - y
      }
      temp <- x; x <- y; y <- temp
    }
    s <- s / 2
  }
  return(d)
}

getHilbertOrder <- function(n) {
  m <- ceiling(log2(n))
  grid_size <- 2^m
  
  coords <- expand.grid(x = 0:(n - 1), y = 0:(n - 1))
  
  coords$x_scaled <- floor(coords$x / n * grid_size)
  coords$y_scaled <- floor(coords$y / n * grid_size)
  
  coords$h_index <- mapply(function(x, y) xy2d(grid_size, x, y), 
                           coords$x_scaled, coords$y_scaled)
  
  order_index <- order(coords$h_index)
  return(order_index)
}

hilbert_order <- getHilbertOrder(n)

stretch_image <- function(image_matrix, order_vec) {
  vec <- as.vector(t(image_matrix))
  
  stretched <- vec[order_vec]
  return(stretched)
}

vessel_hilbert_lines <- vector("list", 3 * N_img_per_group)
maldi_hilbert_lines <- vector("list", 3 * N_img_per_group)

for(i in 1:(3 * N_img_per_group)) {
  vessel_hilbert_lines[[i]] <- stretch_image(vessel_images[[i]], hilbert_order)
  maldi_hilbert_lines[[i]] <- stretch_image(maldi_images[[i]], hilbert_order)
}


reg_data <- data.frame(
  vessel = unlist(vessel_hilbert_lines),
  maldi = unlist(maldi_hilbert_lines),
  group = rep(rep(c("Strong", "Weak", "None"), each = N_img_per_group), each = n*n),
  image_id = rep(1:(3 * N_img_per_group), each = n*n)
)

lm_strong <- lm(maldi ~ vessel, data = subset(reg_data, group == "Strong"))
lm_weak <- lm(maldi ~ vessel, data = subset(reg_data, group == "Weak"))
lm_none <- lm(maldi ~ vessel, data = subset(reg_data, group == "None"))

lm_interaction <- lm(maldi ~ vessel * group, data = reg_data)

cat("\n===== Group-Specific Linear Regression Models =====\n")
cat("\n----- Strong Correlation Group -----\n")
print(summary(lm_strong))

cat("\n----- Weak Correlation Group -----\n")
print(summary(lm_weak))

cat("\n----- No Correlation Group -----\n")
print(summary(lm_none))

cat("\n===== Interaction Model to Test Group Differences =====\n")
print(summary(lm_interaction))
print(anova(lm_interaction))


# 5. ROC Analysis for Discrimination


# Function to calculate regression strength for each image
calculate_r_squared <- function() {
  r_squared_values <- numeric(3 * N_img_per_group)
  p_values <- numeric(3 * N_img_per_group)
  coefficients <- numeric(3 * N_img_per_group)
  
  for(i in 1:(3 * N_img_per_group)) {
    img_data <- data.frame(
      vessel = vessel_hilbert_lines[[i]],
      maldi = maldi_hilbert_lines[[i]]
    )
    lm_img <- lm(maldi ~ vessel, data = img_data)
    summary_lm <- summary(lm_img)
    
    r_squared_values[i] <- summary_lm$r.squared
    p_values[i] <- ifelse(length(summary_lm$coefficients[,4]) > 1, 
                          summary_lm$coefficients[2,4], 1)  
    coefficients[i] <- ifelse(length(summary_lm$coefficients[,1]) > 1,
                              summary_lm$coefficients[2,1], 0) 
  }
  
  return(data.frame(
    image_id = 1:(3 * N_img_per_group),
    group = image_groups,
    group_name = c("Strong", "Weak", "None")[image_groups],
    r_squared = r_squared_values,
    p_value = p_values,
    coefficient = coefficients
  ))
}

image_metrics <- calculate_r_squared()

cat("\n===== Image-Level Regression Metrics =====\n")
print(image_metrics)

group_metrics <- aggregate(
  cbind(r_squared, p_value, coefficient) ~ group_name, 
  data = image_metrics, 
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)

cat("\n===== Group-Level Metrics (Mean ± SD) =====\n")
print(group_metrics)

roc_strong_weak <- roc(
  response = c(rep(1, N_img_per_group), rep(0, N_img_per_group)),
  predictor = image_metrics$r_squared[image_metrics$group %in% c(1, 2)]
)

# Strong vs None
roc_strong_none <- roc(
  response = c(rep(1, N_img_per_group), rep(0, N_img_per_group)),
  predictor = image_metrics$r_squared[image_metrics$group %in% c(1, 3)]
)

# Weak vs None
roc_weak_none <- roc(
  response = c(rep(1, N_img_per_group), rep(0, N_img_per_group)),
  predictor = image_metrics$r_squared[image_metrics$group %in% c(2, 3)]
)

cat("\n===== ROC Analysis for Group Discrimination =====\n")
cat("\nStrong vs Weak: AUC =", auc(roc_strong_weak), "\n")
cat("Strong vs None: AUC =", auc(roc_strong_none), "\n")
cat("Weak vs None: AUC =", auc(roc_weak_none), "\n")



p_rsq <- ggplot(image_metrics, aes(x = factor(group_name, levels = c("Strong", "Weak", "None")), 
                                   y = r_squared, color = group_name)) +
  geom_boxplot(alpha = 0.3) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "R-squared by Group", x = "Group", y = "R-squared") +
  theme_minimal()

p_coef <- ggplot(image_metrics, aes(x = factor(group_name, levels = c("Strong", "Weak", "None")), 
                                    y = coefficient, color = group_name)) +
  geom_boxplot(alpha = 0.3) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "Regression Coefficients by Group", x = "Group", y = "Coefficient") +
  theme_minimal()

plot_regression_example <- function(group_idx = 1) {
  group_indices <- which(image_groups == group_idx)[1]
  
  plot_data <- data.frame(
    vessel = vessel_hilbert_lines[[group_indices]],
    maldi = maldi_hilbert_lines[[group_indices]]
  )
  
  lm_model <- lm(maldi ~ vessel, data = plot_data)
  
  p <- ggplot(plot_data, aes(x = vessel, y = maldi)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, color = "red") +
    labs(title = paste("Regression Example -", 
                       c("Strong Correlation", "Weak Correlation", "No Correlation")[group_idx]),
         subtitle = paste("R² =", round(summary(lm_model)$r.squared, 3),
                          "| p-value =", ifelse(length(summary(lm_model)$coefficients[,4]) > 1,
                                                format.pval(summary(lm_model)$coefficients[2,4], digits = 3), "NA")),
         x = "Vessel (Hilbert-ordered)", y = "MALDI (Hilbert-ordered)") +
    theme_minimal()
  
  return(p)
}

p_strong <- plot_regression_example(1)
p_weak <- plot_regression_example(2)
p_none <- plot_regression_example(3)

p_roc <- ggplot() +
  geom_line(aes(x = roc_strong_weak$specificities, y = roc_strong_weak$sensitivities, 
                color = "Strong vs Weak")) +
  geom_line(aes(x = roc_strong_none$specificities, y = roc_strong_none$sensitivities, 
                color = "Strong vs None")) +
  geom_line(aes(x = roc_weak_none$specificities, y = roc_weak_none$sensitivities, 
                color = "Weak vs None")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves for Group Discrimination", 
       subtitle = paste("AUC: Strong vs Weak =", round(auc(roc_strong_weak), 3),
                        "| Strong vs None =", round(auc(roc_strong_none), 3),
                        "| Weak vs None =", round(auc(roc_weak_none), 3)),
       x = "Specificity", y = "Sensitivity", color = "Comparison") +
  scale_x_reverse() +  
  theme_minimal()

print(p_rsq)
print(p_coef)
print(p_strong)
print(p_weak)
print(p_none)
print(p_roc)



discrimination_assessment <- ifelse(
  min(auc(roc_strong_weak), auc(roc_strong_none), auc(roc_weak_none)) > 0.8,
  "EXCELLENT",
  ifelse(
    min(auc(roc_strong_weak), auc(roc_strong_none), auc(roc_weak_none)) > 0.7,
    "GOOD",
    "MODERATE"
  )
)

cat("Overall discriminative ability of Hilbert Regression method:", discrimination_assessment, "\n")
