# Load the rgl package
library(rgl)

# Clear any previous plots
clear3d()

# Recursive function to draw 3D branches
draw_branch_3D <- function(x, y, z, length, angle_x, angle_y, depth, shrink = 0.7, angle_split = 30) {
  # Base case: stop if depth is zero
  if (depth == 0) return()
  
  # Calculate the end point of the current branch
  x_end <- x + length * sin(angle_y * pi / 180) * cos(angle_x * pi / 180)
  y_end <- y + length * sin(angle_y * pi / 180) * sin(angle_x * pi / 180)
  z_end <- z + length * cos(angle_y * pi / 180)
  
  # Draw the branch
  segments3d(rbind(c(x, y, z), c(x_end, y_end, z_end)), col = "brown", lwd = depth)
  
  # Recursive calls for the branches
  # Left branch
  draw_branch_3D(x_end, y_end, z_end, length * shrink, angle_x - angle_split, angle_y - angle_split, depth - 1, shrink, angle_split)
  
  # Right branch
  draw_branch_3D(x_end, y_end, z_end, length * shrink, angle_x + angle_split, angle_y - angle_split, depth - 1, shrink, angle_split)
  
  # Upward branch
  draw_branch_3D(x_end, y_end, z_end, length * shrink, angle_x, angle_y + angle_split, depth - 1, shrink, angle_split)
}

# Initialize 3D plot
open3d()
bg3d("white")  # Set background color
view3d(theta = 30, phi = 20, zoom = 0.7)  # Set viewing angle

# Set initial parameters for the tree structure
trunk_length <- 5
initial_x <- 0
initial_y <- 0
initial_z <- 0
initial_angle_x <- 0   # Horizontal angle
initial_angle_y <- 90  # Vertical angle
depth <- 6             # Recursion depth

# Draw the 3D fractal tree structure
draw_branch_3D(initial_x, initial_y, initial_z, trunk_length, initial_angle_x, initial_angle_y, depth)

# Add axes for reference
axes3d()
title3d("3D Fractal vessel", "", "", "", "")
##################################################
# Clear any previous plots
clear3d()


# Define a 3D space as boundaries
limi <- 100
space_x <- c(-limi, limi)  # X bounds
space_y <- c(-limi, limi)  # Y bounds
space_z <- c(-limi, limi)    # Z bounds

# Recursive function to simulate vascular-like branches
draw_vasculature_dense <- function(x, y, z, length, angle_x, angle_y, depth, 
                                   shrink = 0.75, angle_variation = 25, space_limits) {
  # Base case: stop if depth is zero or length is very small
  if (depth == 0 || length < 0.3) return()
  
  # Calculate the end point of the current branch
  x_end <- x + length * sin(angle_y * pi / 180) * cos(angle_x * pi / 180)
  y_end <- y + length * sin(angle_y * pi / 180) * sin(angle_x * pi / 180)
  z_end <- z + length * cos(angle_y * pi / 180)
  
  # Check if the branch end is outside the defined space
  if (x_end < space_limits$x[1] || x_end > space_limits$x[2] || 
      y_end < space_limits$y[1] || y_end > space_limits$y[2] || 
      z_end < space_limits$z[1] || z_end > space_limits$z[2]) {
    return()
  }
  
  # Draw the current branch
  segments3d(rbind(c(x, y, z), c(x_end, y_end, z_end)), col = "red", lwd = max(1, depth / 2))
  
  # Introduce randomness in the angles for a natural look
  for (i in 1:4) {  # Increase the number of branches per node
    new_angle_x <- angle_x + runif(1, -angle_variation, angle_variation)
    new_angle_y <- angle_y + runif(1, -angle_variation, angle_variation)
    draw_vasculature_dense(x_end, y_end, z_end, length * shrink, new_angle_x, new_angle_y, 
                           depth - 1, shrink, angle_variation, space_limits)
  }
}

# Define the space boundaries
space_limits <- list(x = space_x, y = space_y, z = space_z)

# Initialize the 3D plot
open3d()
bg3d("white")  # Set background to white
view3d(theta = 30, phi = 40, zoom = 0.7)  # Adjust the viewpoint

# Draw the initial trunk/vasculature root
initial_length <- 3
initial_depth <- 7  # Increased recursion depth
# Add axes and title
axes3d()
title3d("Dense 3D Tissue Vasculature Simulation", "", "", "", "")
# Start recursive function at the base of the space
draw_vasculature_dense(0, 0, 0, initial_length, angle_x = 0, angle_y = 90, 
                       depth = initial_depth, space_limits = space_limits)
# 
# # Draw a bounding box to represent the tissue space
# lines3d(rbind(
#   c(space_x[1], space_y[1], space_z[1]),
#   c(space_x[2], space_y[1], space_z[1]),
#   c(space_x[2], space_y[2], space_z[1]),
#   c(space_x[1], space_y[2], space_z[1]),
#   c(space_x[1], space_y[1], space_z[1]),
#   c(space_x[1], space_y[1], space_z[2]),
#   c(space_x[2], space_y[1], space_z[2]),
#   c(space_x[2], space_y[2], space_z[2]),
#   c(space_x[1], space_y[2], space_z[2]),
#   c(space_x[1], space_y[1], space_z[2])
# ), col = "black", lwd = 2)

# Save the 3D object as an .obj file
writeOBJ("output/vasculature.obj")

rgl.write("vasculature.obj")

# Load a 3D mesh object (e.g., from an .obj file)
mesh <- readOBJ("output/vasculature.obj")

# Extract vertices (3D coordinates of the object)
vertices <- mesh$vb

# Define the slicing plane (along the Z-axis)
slice_plane_z <- 0
tolerance <- 0.1  # Range around the plane for slicing

# Extract the points near the slicing plane
slice_points <- vertices[, abs(vertices[3, ] - slice_plane_z) < tolerance]

# Plot the 2D slice (project the X and Y coordinates)
plot(slice_points[1, ], slice_points[2, ], pch = 20, col = "blue", 
     main = "2D Slice of 3D Mesh", xlab = "X", ylab = "Y")
