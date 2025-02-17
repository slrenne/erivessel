# Load necessary library
set.seed(123)  # For reproducibility

# Step 1: Generate multiple random polygons
generate_random_polygon <- function(center_x, center_y, num_sides, radius) {
  angles <- seq(0, 2 * pi, length.out = num_sides + 1)
  x <- center_x + radius * cos(angles)
  y <- center_y + radius * sin(angles)
  return(data.frame(x = x, y = y))
}

# Create random polygons
polygon1 <- generate_random_polygon(2, 2, 4, 1)  # Square near (2, 2)
polygon2 <- generate_random_polygon(5, 4, 5, 1.5)  # Pentagon near (5, 4)
polygon3 <- generate_random_polygon(7, 6, 6, 1)  # Hexagon near (7, 6)
polygon4 <- generate_random_polygon(3, 5, 4, 0.8)  # Smaller square near (3, 5)

# Step 2: Plot the polygons
plot(0, 0, type = "n", xlim = c(0, 8), ylim = c(0, 8), xlab = "X-Axis", ylab = "Y-Axis",
     main = "Polygons and Multiple Random Lines")

polygon(polygon1$x, polygon1$y, col = "lightblue", border = "blue")
polygon(polygon2$x, polygon2$y, col = "lightgreen", border = "darkgreen")
polygon(polygon3$x, polygon3$y, col = "lightpink", border = "red")
polygon(polygon4$x, polygon4$y, col = "lightyellow", border = "orange")

# Step 3: Function to check if two line segments intersect
do_segments_intersect <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
  # Calculate the direction of each segment
  det <- (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3)
  if (det == 0) return(FALSE)  # Parallel or collinear segments
  
  # Calculate intersection point using parametric equations
  t <- ((x3 - x1) * (y4 - y3) - (y3 - y1) * (x4 - x3)) / det
  u <- ((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)) / det
  
  # Check if intersection lies on both segments
  return(t >= 0 && t <= 1 && u >= 0 && u <= 1)
}

# Step 4: Function to draw a random line and highlight overlapping segments
draw_random_line_and_highlight_overlaps <- function() {
  slope <- runif(1, -2, 2)   # Random slope between -2 and 2
  intercept <- runif(1, -5, 5)  # Random intercept between -5 and 5
  line_x <- seq(0, 8, by = 0.01)
  line_y <- slope * line_x + intercept
  
  # Draw the random line
  lines(line_x, line_y, col = "red", lwd = 2)
  
  # Step 5: Highlight overlapping segments
  for (polygon in list(polygon1, polygon2, polygon3, polygon4)) {
    for (i in 1:(length(line_x) - 1)) {
      for (j in 1:(nrow(polygon) - 1)) {
        if (do_segments_intersect(line_x[i], line_y[i], line_x[i + 1], line_y[i + 1],
                                  polygon$x[j], polygon$y[j], polygon$x[j + 1], polygon$y[j + 1])) {
          # Highlight the entire overlapping segment
          segments(line_x[i], line_y[i], line_x[i + 1], line_y[i + 1], col = "blue", lwd = 3)
        }
      }
    }
  }
}

# Step 6: Draw multiple random lines and highlight overlaps
num_lines <- 5  # Number of random lines to draw
for (i in 1:num_lines) {
  draw_random_line_and_highlight_overlaps()
}
