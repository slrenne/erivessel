library(ggplot2)
library(dplyr)
library(FSA)      # For Dunn's test
library(effsize)  # For effect size calculations
library(nortest)  # For Anderson-Darling normality test

# Read data
my_data <- read.csv("db.csv", stringsAsFactors = FALSE)

# Create directories for output
if(!dir.exists("plots")) dir.create("plots")
if(!dir.exists("results")) dir.create("results")

# Get unique models
models <- unique(my_data$model)

# Create results file
results_file <- file("results/statistical_results.txt", "w")
cat("Statistical Analysis Results\n\n", file = results_file)

# Function to handle large sample normality testing
test_normality <- function(x) {
  if(length(x) > 5000) {
    list(
      test = "Anderson-Darling (n > 5000)",
      p.value = ad.test(x)$p.value
    )
  } else {
    list(
      test = "Shapiro-Wilk",
      p.value = shapiro.test(x)$p.value
    )
  }
}

for (i in seq_along(models)) {
  current_model <- models[i]
  model_data <- my_data %>% filter(model == current_model)
  
  # Create enhanced plot
  p <- ggplot(model_data, aes(x = axis_minor_length, color = treatment)) +
    geom_freqpoly(bins = 30, linewidth = 1) +  # Changed from size to linewidth
    geom_boxplot(aes(y = 0, fill = treatment), width = 2, alpha = 0.3, 
                 position = position_dodge(width = 0.5)) +
    stat_summary(aes(y = 0), fun.data = "mean_cl_normal", 
                 geom = "pointrange", color = "black", linewidth = 0.5,
                 na.rm = TRUE) +  # Added na.rm
    labs(title = paste("Model:", current_model),
         x = "Axis Minor Length",
         y = "Number of Vessels",
         color = "Treatment",
         fill = "Treatment") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Display and save plot
  print(p)
  ggsave(
    filename = paste0("plots/model_", current_model, ".png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # Statistical analysis
  cat("\n\n===========================================\n", 
      file = results_file, append = TRUE)
  cat(paste("Statistical Analysis for Model:", current_model, "\n"), 
      file = results_file, append = TRUE)
  cat("===========================================\n\n", 
      file = results_file, append = TRUE)
  
  treatments <- unique(model_data$treatment)
  num_treatments <- length(treatments)
  
  if(num_treatments > 1) {
    # Descriptive statistics
    desc_stats <- model_data %>%
      group_by(treatment) %>%
      summarise(
        n = n(),
        mean = mean(axis_minor_length, na.rm = TRUE),
        sd = sd(axis_minor_length, na.rm = TRUE),
        median = median(axis_minor_length, na.rm = TRUE),
        IQR = IQR(axis_minor_length, na.rm = TRUE)
      )
    
    cat("Descriptive Statistics:\n", file = results_file, append = TRUE)
    capture.output(print(desc_stats), file = results_file, append = TRUE)
    cat("\n", file = results_file, append = TRUE)
    
    # Normality assessment
    cat("Normality Assessment:\n", file = results_file, append = TRUE)
    norm_results <- by(model_data$axis_minor_length, model_data$treatment, test_normality)
    capture.output(print(norm_results), file = results_file, append = TRUE)
    cat("\n", file = results_file, append = TRUE)
    
    # Determine if all groups are normal
    all_normal <- all(sapply(norm_results, function(x) x$p.value > 0.05))
    
    if(num_treatments > 2) {
      # Multiple groups comparison
      if(all_normal) {
        # ANOVA
        anova_result <- aov(axis_minor_length ~ treatment, data = model_data)
        cat("ANOVA Results:\n", file = results_file, append = TRUE)
        capture.output(print(summary(anova_result)), file = results_file, append = TRUE)
        
        if(summary(anova_result)[[1]]$`Pr(>F)`[1] < 0.05) {
          # Tukey HSD
          cat("\nTukey HSD Post-hoc Test:\n", file = results_file, append = TRUE)
          tukey_result <- TukeyHSD(anova_result)
          capture.output(print(tukey_result), file = results_file, append = TRUE)
        }
      } else {
        # Kruskal-Wallis
        kruskal_result <- kruskal.test(axis_minor_length ~ treatment, data = model_data)
        cat("Kruskal-Wallis Test Results:\n", file = results_file, append = TRUE)
        capture.output(print(kruskal_result), file = results_file, append = TRUE)
        
        if(kruskal_result$p.value < 0.05) {
          # Dunn's test
          cat("\nDunn's Post-hoc Test:\n", file = results_file, append = TRUE)
          dunn_result <- dunnTest(axis_minor_length ~ treatment, 
                                  data = model_data, method = "bh")
          capture.output(print(dunn_result), file = results_file, append = TRUE)
        }
      }
    } else {
      # Two groups comparison
      group1 <- model_data$axis_minor_length[model_data$treatment == treatments[1]]
      group2 <- model_data$axis_minor_length[model_data$treatment == treatments[2]]
      
      if(all_normal) {
        # t-test
        t_test_result <- t.test(group1, group2)
        cat("Student's t-test Results:\n", file = results_file, append = TRUE)
        capture.output(print(t_test_result), file = results_file, append = TRUE)
      } else {
        # Wilcoxon test
        wilcox_result <- wilcox.test(group1, group2)
        cat("Wilcoxon Rank Sum Test Results:\n", file = results_file, append = TRUE)
        capture.output(print(wilcox_result), file = results_file, append = TRUE)
      }
    }
    
    # Effect sizes
    cat("\nEffect Sizes:\n", file = results_file, append = TRUE)
    treatment_pairs <- combn(treatments, 2, simplify = FALSE)
    
    for(pair in treatment_pairs) {
      g1 <- model_data$axis_minor_length[model_data$treatment == pair[1]]
      g2 <- model_data$axis_minor_length[model_data$treatment == pair[2]]
      
      if(all_normal) {
        es <- cohen.d(g1, g2)
        cat(paste("\nCohen's d for", pair[1], "vs", pair[2], ":\n"), 
            file = results_file, append = TRUE)
      } else {
        es <- cliff.delta(g1, g2)
        cat(paste("\nCliff's delta for", pair[1], "vs", pair[2], ":\n"), 
            file = results_file, append = TRUE)
      }
      capture.output(print(es), file = results_file, append = TRUE)
    }
  } else {
    cat("Only one treatment group found - no comparisons to make.\n", 
        file = results_file, append = TRUE)
  }
  
  cat("\nProcessed model:", current_model, "\n")
}

close(results_file)

cat("\nAnalysis complete!\n")
cat("- Plots saved in 'plots/' directory\n")
cat("- Statistical results saved in 'results/statistical_results.txt'\n")
