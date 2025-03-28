library(ggplot2)
library(dplyr)

# Read data
data <- read.csv("dbsim.csv", stringsAsFactors = FALSE)

# Calculate summary statistics
stats <- data %>%
  summarise(
    Mean = mean(minor_axis),
    Median = median(minor_axis),
    SD = sd(minor_axis),
    Min = min(minor_axis),
    Max = max(minor_axis),
    N = n()
  ) %>%
  mutate(Label = sprintf("Mean: %.3f\nMedian: %.3f\nSD: %.3f\nMin: %.3f\nMax: %.3f\nN: %d",
                         Mean, Median, SD, Min, Max, N))

# Create plots directory if it doesn't exist
if(!dir.exists("plots")) {
  dir.create("plots")
}

# Create the frequency plot with statistics
p <- ggplot(data, aes(x = minor_axis)) +
  geom_histogram(bins = 30, fill = "#56B4E9", alpha = 0.8, color = "white") +
  geom_freqpoly(bins = 30, linewidth = 1, color = "#0072B2") +
  geom_vline(xintercept = stats$Mean, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  geom_vline(xintercept = stats$Median, linetype = "dashed", color = "#009E73", linewidth = 0.8) +
  annotate("text", 
           x = Inf, y = Inf,
           label = stats$Label,
           hjust = 1.1, vjust = 1.1,
           size = 4, color = "black") +
  labs(title = "Distribution of Minor Axis Lengths",
       x = "Minor Axis Length",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12))

# Display plot
print(p)

# Save plot
ggsave(
  filename = "plots/minor_axis_distribution_with_stats.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Plot saved as 'plots/minor_axis_distribution_with_stats.png'\n")