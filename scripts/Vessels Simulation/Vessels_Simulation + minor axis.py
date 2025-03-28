import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse
import random
from dataclasses import dataclass
import math
import time
from matplotlib.animation import FuncAnimation
from IPython.display import clear_output
import csv  # Added for CSV export
import os  # Added for file path handling

@dataclass
class Vessel:
    start: np.ndarray  # 3D point where vessel starts
    end: np.ndarray    # 3D point where vessel ends
    radius: float      # vessel radius
    parent: object = None  # parent vessel

class VascularSystem:
    def __init__(self, start_point, direction, initial_length, initial_radius, 
                 min_branching_angle=10, max_branching_angle=70, radius_ratio=0.5, length_ratio=0.6, 
                 max_depth=9, min_radius=0.1, random_seed=None, 
                 progressive_growth=True, visualization_interval=30.0):
        """
        Initialize vascular system with fractal-like bifurcations
        
        Parameters:
        - start_point: 3D starting point of the root vessel
        - direction: Initial direction of growth
        - initial_length: Length of the root vessel
        - initial_radius: Radius of the root vessel
        - min_branching_angle: Minimum branching angle in degrees
        - max_branching_angle: Maximum branching angle in degrees
        - radius_ratio: Ratio of child radius to parent radius
        - length_ratio: Ratio of child length to parent length
        - max_depth: Maximum recursion depth for bifurcations
        - min_radius: Minimum vessel radius to stop bifurcation
        - random_seed: Seed for random number generation (for reproducibility)
        - progressive_growth: Whether to visualize growth progressively
        - visualization_interval: Time interval (in seconds) between progressive visualizations
        """
        # Set random seed if provided
        if random_seed is not None:
            random.seed(random_seed)
            np.random.seed(random_seed)
            
        self.vessels = []
        self.start_point = np.array(start_point)
        self.direction = self._normalize(np.array(direction))
        self.initial_length = initial_length
        self.initial_radius = initial_radius
        self.min_branching_angle = np.radians(min_branching_angle)
        self.max_branching_angle = np.radians(max_branching_angle)
        self.radius_ratio = radius_ratio
        self.length_ratio = length_ratio
        self.max_depth = max_depth
        self.min_radius = min_radius
        self.progressive_growth = progressive_growth
        self.visualization_interval = visualization_interval
        
        # For progressive growth visualization
        self.growth_fig = None
        self.growth_ax = None
        if progressive_growth:
            self.growth_fig = plt.figure(figsize=(10, 8))
            self.growth_ax = self.growth_fig.add_subplot(111, projection='3d')
            
        # Create the root vessel
        self._create_root_vessel()
        
        # Vessels organized by generation for progressive growth
        self.vessels_by_generation = {0: [self.vessels[0]]}
        
        # Generate the fractal vascular system
        if progressive_growth:
            self._generate_vessels_progressively(self.vessels[0], 1)
        else:
            self._generate_vessels(self.vessels[0], 1)
    
    def _normalize(self, vector):
        """Normalize a vector to unit length"""
        norm = np.linalg.norm(vector)
        if norm == 0:
            return vector
        return vector / norm
    
    def _create_root_vessel(self):
        """Create the initial root vessel"""
        end_point = self.start_point + self.direction * self.initial_length
        root_vessel = Vessel(
            start=self.start_point,
            end=end_point,
            radius=self.initial_radius
        )
        self.vessels.append(root_vessel)
    
    def _get_perpendicular_vector(self, vector):
        """Find a perpendicular vector to the given vector"""
        # Create a non-parallel vector
        if abs(vector[0]) < abs(vector[1]):
            perpendicular = np.array([1, 0, 0])
        else:
            perpendicular = np.array([0, 1, 0])
            
        # Use cross product to get perpendicular vector
        return self._normalize(np.cross(vector, perpendicular))
    
    def _rotate_vector(self, vector, axis, angle):
        """Rotate a vector around an axis by the given angle"""
        # Rodrigues rotation formula
        axis = self._normalize(axis)
        k_cross_v = np.cross(axis, vector)
        return (vector * np.cos(angle) + 
                k_cross_v * np.sin(angle) + 
                axis * np.dot(axis, vector) * (1 - np.cos(angle)))
    
    def _generate_vessels(self, parent_vessel, depth):
        """Recursively generate vessels with bifurcations"""
        if depth >= self.max_depth or parent_vessel.radius <= self.min_radius:
            return
        
        # Calculate child vessel parameters
        parent_direction = self._normalize(parent_vessel.end - parent_vessel.start)
        child_radius = parent_vessel.radius * self.radius_ratio
        child_length = self.initial_length * (self.length_ratio ** depth)
        
        # Get a perpendicular axis for rotation
        rotation_axis = self._get_perpendicular_vector(parent_direction)
        
        # Increase branching factor as we go deeper (more bifurcations for smaller vessels)
        branching_factor = 2
        if depth > 2:  # After certain depth, increase branching
            # More branches at deeper levels, max 4 branches
            branching_factor = min(4, 2 + int(depth/2))
        
        # Distribute branches around the parent vessel
        angle_step = 2 * np.pi / branching_factor
        
        # Create branches at angles from the parent direction
        for i in range(branching_factor):
            # Distribute branches evenly around the vessel
            azimuth_angle = i * angle_step
            
            # Add some random variation to the azimuth angle
            azimuth_variation = np.radians(random.uniform(-15, 15))
            actual_azimuth = azimuth_angle + azimuth_variation
            
            # Randomize the branching angle within the specified range
            actual_angle = random.uniform(self.min_branching_angle, self.max_branching_angle)
            
            # First rotate around parent direction to set azimuth
            azimuth_rotation_axis = parent_direction
            temp_axis = self._rotate_vector(rotation_axis, azimuth_rotation_axis, actual_azimuth)
            
            # Then rotate around the temporary axis for branching angle
            child_direction = self._rotate_vector(parent_direction, temp_axis, actual_angle)
            
            # Calculate child vessel endpoints
            child_start = parent_vessel.end
            child_end = child_start + child_direction * child_length
            
            # Create the child vessel
            child_vessel = Vessel(
                start=child_start,
                end=child_end,
                radius=child_radius,
                parent=parent_vessel
            )
            
            # Add the vessel to our collection
            self.vessels.append(child_vessel)
            
            # Add to vessels by generation dictionary
            if depth not in self.vessels_by_generation:
                self.vessels_by_generation[depth] = []
            self.vessels_by_generation[depth].append(child_vessel)
            
            # Recursively generate the next level of vessels
            self._generate_vessels(child_vessel, depth + 1)
            
    def _generate_vessels_progressively(self, parent_vessel, depth):
        """Generate vessels with progressive visualization"""
        if depth >= self.max_depth or parent_vessel.radius <= self.min_radius:
            return
        
        # Generate vessels for this generation
        current_generation_vessels = []
        
        # Calculate child vessel parameters
        parent_direction = self._normalize(parent_vessel.end - parent_vessel.start)
        child_radius = parent_vessel.radius * self.radius_ratio
        child_length = self.initial_length * (self.length_ratio ** depth)
        
        # Get a perpendicular axis for rotation
        rotation_axis = self._get_perpendicular_vector(parent_direction)
        
        # Increase branching factor as we go deeper
        branching_factor = 2
        if depth > 2:  # After certain depth, increase branching
            branching_factor = min(4, 2 + int(depth/2))
        
        # Distribute branches around the parent vessel
        angle_step = 2 * np.pi / branching_factor
        
        # Create branches at angles from the parent direction
        for i in range(branching_factor):
            azimuth_angle = i * angle_step
            azimuth_variation = np.radians(random.uniform(-15, 15))
            actual_azimuth = azimuth_angle + azimuth_variation
            
            # Randomize the branching angle within the specified range
            actual_angle = random.uniform(self.min_branching_angle, self.max_branching_angle)
            
            azimuth_rotation_axis = parent_direction
            temp_axis = self._rotate_vector(rotation_axis, azimuth_rotation_axis, actual_azimuth)
            
            child_direction = self._rotate_vector(parent_direction, temp_axis, actual_angle)
            
            child_start = parent_vessel.end
            child_end = child_start + child_direction * child_length
            
            child_vessel = Vessel(
                start=child_start,
                end=child_end,
                radius=child_radius,
                parent=parent_vessel
            )
            
            self.vessels.append(child_vessel)
            current_generation_vessels.append(child_vessel)
            
            # Add to vessels by generation dictionary
            if depth not in self.vessels_by_generation:
                self.vessels_by_generation[depth] = []
            self.vessels_by_generation[depth].append(child_vessel)
        
        # Visualize this generation
        if self.progressive_growth and self.growth_ax:
            self._visualize_generation(depth)
            plt.pause(self.visualization_interval)
        
        # Recursively generate the next level for each child
        for child_vessel in current_generation_vessels:
            self._generate_vessels_progressively(child_vessel, depth + 1)
            
    def _visualize_generation(self, current_depth):
        """Visualize vessels up to the current generation"""
        self.growth_ax.clear()
        
        # Plot vessels up to the current generation
        for depth in range(current_depth + 1):
            if depth in self.vessels_by_generation:
                for vessel in self.vessels_by_generation[depth]:
                    # Use different colors for different generations
                    color = plt.cm.viridis(depth / self.max_depth)
                    linewidth = max(0.5, vessel.radius * 3)  # Scale linewidth
                    
                    self.growth_ax.plot(
                        [vessel.start[0], vessel.end[0]],
                        [vessel.start[1], vessel.end[1]],
                        [vessel.start[2], vessel.end[2]],
                        color=color, linewidth=linewidth
                    )
        
        # Set labels and title
        self.growth_ax.set_xlabel('X')
        self.growth_ax.set_ylabel('Y')
        self.growth_ax.set_zlabel('Z')
        self.growth_ax.set_title(f'Vascular System Growth - Generation {current_depth}')
        
        # Auto-adjust limits with 10% margin
        all_points = []
        for vessels in self.vessels_by_generation.values():
            for v in vessels:
                all_points.append(v.start)
                all_points.append(v.end)
        
        if all_points:
            all_points = np.array(all_points)
            min_vals = np.min(all_points, axis=0)
            max_vals = np.max(all_points, axis=0)
            ranges = max_vals - min_vals
            margin = 0.1 * np.max(ranges)
            
            self.growth_ax.set_xlim(min_vals[0] - margin, max_vals[0] + margin)
            self.growth_ax.set_ylim(min_vals[1] - margin, max_vals[1] + margin)
            self.growth_ax.set_zlim(min_vals[2] - margin, max_vals[2] + margin)
        
        plt.draw()
    
    def visualize_3d(self, ax=None, show_plane=False, plane_params=None):
        """Visualize the vascular system in 3D"""
        if ax is None:
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(111, projection='3d')
            
        # Plot each vessel as a line segment
        for vessel in self.vessels:
            ax.plot(
                [vessel.start[0], vessel.end[0]],
                [vessel.start[1], vessel.end[1]],
                [vessel.start[2], vessel.end[2]],
                'r-', linewidth=vessel.radius * 5  # Scale the linewidth
            )
        
        # If provided, show the intersection plane
        if show_plane and plane_params:
            self._plot_plane(ax, plane_params)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D Vascular System')
        
        return ax
    
    def _plot_plane(self, ax, plane_params):
        """Plot the intersection plane in 3D"""
        a, b, c, d = plane_params  # ax + by + cz + d = 0
        
        # Find the limits of our vascular system
        xs = [v.start[0] for v in self.vessels] + [v.end[0] for v in self.vessels]
        ys = [v.start[1] for v in self.vessels] + [v.end[1] for v in self.vessels]
        zs = [v.start[2] for v in self.vessels] + [v.end[2] for v in self.vessels]
        
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        margin = 0.2 * max(max_x - min_x, max_y - min_y)
        
        # Extend boundaries with margins
        x_range = np.linspace(min_x - margin, max_x + margin, 10)
        y_range = np.linspace(min_y - margin, max_y + margin, 10)
        xx, yy = np.meshgrid(x_range, y_range)
        
        # Calculate z from the plane equation: ax + by + cz + d = 0
        if c != 0:  # Make sure c is not zero to avoid division by zero
            zz = (-a * xx - b * yy - d) / c
            
            # Plot the plane
            ax.plot_surface(xx, yy, zz, alpha=0.3, color='blue')
        else:
            # Handle case where plane is parallel to z-axis (c=0)
            min_z, max_z = min(zs), max(zs)
            z_range = np.linspace(min_z - margin, max_z + margin, 10)
            
            if a != 0:  # Plane is of form ax + by + d = 0 (c=0)
                yy, zz = np.meshgrid(y_range, z_range)
                xx = (-b * yy - d) / a
                ax.plot_surface(xx, yy, zz, alpha=0.3, color='blue')
            elif b != 0:  # Plane is of form by + d = 0 (a=0, c=0)
                xx, zz = np.meshgrid(x_range, z_range)
                yy = -d / b
                ax.plot_surface(xx, yy, zz, alpha=0.3, color='blue')
    
    def generate_random_plane(self):
        """Generate a random plane that intersects with the vascular system"""
        # Find center of the vascular system
        all_points = np.array([v.start for v in self.vessels] + [v.end for v in self.vessels])
        center = np.mean(all_points, axis=0)
        
        # Generate random normal vector for the plane
        normal = np.array([
            random.uniform(-1, 1),
            random.uniform(-1, 1),
            random.uniform(-1, 1)
        ])
        normal = self._normalize(normal)
        
        # Plane equation: ax + by + cz + d = 0
        # Where (a, b, c) is the normal vector
        a, b, c = normal
        # Calculate d so that the plane passes through the center
        d = -(a * center[0] + b * center[1] + c * center[2])
        
        return (a, b, c, d)
    
    def calculate_intersection(self, plane_params):
        """
        Calculate intersections of the plane with all vessels
        
        Returns a list of ellipses (or circles) that represent the intersections
        Each ellipse is described by its center, major and minor axes, and rotation angle
        """
        a, b, c, d = plane_params
        normal = np.array([a, b, c])
        normal = self._normalize(normal)
        
        intersections = []
        
        for vessel in self.vessels:
            # Vessel is represented by a cylinder
            # Get the axis of the cylinder (direction of the vessel)
            vessel_direction = self._normalize(vessel.end - vessel.start)
            
            # Check if the vessel and plane intersect
            p0 = vessel.start
            p1 = vessel.end
            
            # The plane equation: ax + by + cz + d = 0
            # The line equation: p = p0 + t * (p1 - p0), where t is a parameter
            
            # Substituting the line equation into the plane equation:
            # a(p0_x + t(p1_x - p0_x)) + b(p0_y + t(p1_y - p0_y)) + c(p0_z + t(p1_z - p0_z)) + d = 0
            # Solving for t:
            denom = np.dot(normal, p1 - p0)
            
            # If denominator is close to zero, the line is nearly parallel to the plane
            if abs(denom) < 1e-10:
                continue
            
            t = -(np.dot(normal, p0) + d) / denom
            
            # Check if the intersection point is within the vessel segment (0 <= t <= 1)
            if 0 <= t <= 1:
                # Calculate the intersection point
                intersection_point = p0 + t * (p1 - p0)
                
                # The intersection of a plane with a cylinder creates an ellipse
                # The angle between the cylinder axis and the plane normal affects the shape
                cos_angle = abs(np.dot(vessel_direction, normal))
                sin_angle = np.sqrt(1 - cos_angle**2)
                
                # If cylinder axis is perpendicular to plane normal, intersection is a circle
                # Otherwise, it's an ellipse
                if sin_angle < 1e-10:
                    # Circle with radius equal to the vessel radius
                    major_axis = minor_axis = vessel.radius
                    ellipse_angle = 0
                else:
                    # For an ellipse, the major axis is vessel.radius / sin_angle
                    major_axis = vessel.radius / sin_angle
                    minor_axis = vessel.radius
                    
                    # Calculate rotation angle of the ellipse in the plane
                    # Project vessel direction onto the plane
                    proj = vessel_direction - np.dot(vessel_direction, normal) * normal
                    proj = self._normalize(proj)
                    
                    # Define a reference direction in the plane (for angle calculation)
                    ref_dir = self._get_perpendicular_vector(normal)
                    
                    # Calculate the angle between projected direction and reference direction
                    cos_ellipse = np.dot(proj, ref_dir)
                    sin_ellipse = np.dot(np.cross(ref_dir, proj), normal)
                    ellipse_angle = np.arctan2(sin_ellipse, cos_ellipse)
                
                # Store the intersection information
                intersections.append({
                    'center': intersection_point,
                    'major_axis': major_axis,
                    'minor_axis': minor_axis,
                    'angle': ellipse_angle,
                    'normal': normal,
                    'vessel': vessel
                })
        
        return intersections
    
    def save_minor_axes_to_csv(self, intersections, filename="minor_axes.csv"):
        """Save minor axis data to a CSV file"""
        # Get the current working directory
        current_dir = os.getcwd()
        full_path = os.path.join(current_dir, filename)
        
        # Save the file
        with open(full_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['minor_axis'])  # header
            for intersection in intersections:
                writer.writerow([intersection['minor_axis']])
        
        print(f"Minor axis data saved to: {full_path}")
    
    def visualize_intersection(self, plane_params, intersections=None, save_csv=False):
        """Visualize the intersection of the vascular system with a plane"""
        if intersections is None:
            intersections = self.calculate_intersection(plane_params)
        
        if not intersections:
            print("No intersections found!")
            return
        
        # Save minor axes to CSV if requested
        if save_csv:
            self.save_minor_axes_to_csv(intersections)
        
        # Create 3D plot with vessels and plane
        fig = plt.figure(figsize=(16, 8))
        ax1 = fig.add_subplot(121, projection='3d')
        self.visualize_3d(ax1, show_plane=True, plane_params=plane_params)
        
        # Create 2D plot for the intersection
        ax2 = fig.add_subplot(122)
        
        # Project the intersection points onto a 2D coordinate system in the plane
        a, b, c, d = plane_params
        normal = np.array([a, b, c])
        normal = self._normalize(normal)
        
        # Get two orthogonal basis vectors in the plane
        u = self._get_perpendicular_vector(normal)
        v = np.cross(normal, u)
        
        # Origin for the 2D coordinate system (an arbitrary point on the plane)
        if c != 0:
            origin = np.array([0, 0, -d/c])
        elif b != 0:
            origin = np.array([0, -d/b, 0])
        else:
            origin = np.array([-d/a, 0, 0])
        
        # Plot each intersection as an ellipse
        for intersection in intersections:
            center_3d = intersection['center']
            
            # Project the center to 2D coordinates in the plane
            center_vec = center_3d - origin
            center_2d = (np.dot(center_vec, u), np.dot(center_vec, v))
            
            # Create ellipse
            ellipse = Ellipse(
                xy=center_2d,
                width=2 * intersection['major_axis'],
                height=2 * intersection['minor_axis'],
                angle=np.degrees(intersection['angle']),
                edgecolor='red',
                facecolor='none',
                linewidth=2
            )
            
            ax2.add_patch(ellipse)

        
        # Set equal aspect ratio for the 2D plot
        ax2.set_aspect('equal')
        
        # Set limits with some margin
        all_x = [e['center'][0] for e in intersections]
        all_y = [e['center'][1] for e in intersections]
        max_radius = max([max(e['major_axis'], e['minor_axis']) for e in intersections])
        margin = max_radius * 2
        
        min_x, max_x = min(all_x), max(all_x)
        min_y, max_y = min(all_y), max(all_y)
        
        x_range = max_x - min_x
        y_range = max_y - min_y
        
        # Ensure we have some range even if all points are clustered
        if x_range < margin:
            center_x = (min_x + max_x) / 2
            min_x = center_x - margin/2
            max_x = center_x + margin/2
        
        if y_range < margin:
            center_y = (min_y + max_y) / 2
            min_y = center_y - margin/2
            max_y = center_y + margin/2
        
        ax2.set_xlim(min_x - margin, max_x + margin)
        ax2.set_ylim(min_y - margin, max_y + margin)
        
        ax2.set_title('2D Intersection View')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.grid(True)
        
        plt.tight_layout()
        return fig

def create_animation(vascular_system, plane_params=None, save_animation=False):
    """Create an animation of rotating the 3D vascular system"""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Visualize the system
    vascular_system.visualize_3d(ax, show_plane=(plane_params is not None), plane_params=plane_params)
    
    # Set up the animation
    def update(frame):
        ax.view_init(elev=20, azim=frame)
        return fig,
    
    # Create the animation
    ani = FuncAnimation(fig, update, frames=range(0, 360, 2), interval=50, blit=True)
    
    # Save if requested
    if save_animation:
        ani.save('vascular_system_rotation.gif', writer='pillow', fps=20)
    
    plt.tight_layout()
    plt.show()
    
    return ani

def main():
    # Set the random seed for reproducibility
    seed = 333
    
    # Create a vascular system with fractal-like bifurcations
    vascular_system = VascularSystem(
        start_point=[0, 0, 0],
        direction=[0, 0, 1],
        initial_length=10,
        initial_radius=2.0,
        min_branching_angle=10,  # Minimum branching angle
        max_branching_angle=70,  # Maximum branching angle
        radius_ratio=0.5,
        length_ratio=0.6,
        max_depth=9,  # Increased depth for more bifurcations
        min_radius=0.01,  # Smaller minimum radius to allow more bifurcations
        random_seed=seed,
        progressive_growth=True,  # Show growth process
        visualization_interval=0.01  # Time between visualizations
    )
    
    # After growth is complete, generate a random plane
    plane_params = vascular_system.generate_random_plane()
    print(f"Random plane equation: {plane_params[0]}x + {plane_params[1]}y + {plane_params[2]}z + {plane_params[3]} = 0")
    
    # Calculate intersections with the plane
    intersections = vascular_system.calculate_intersection(plane_params)
    print(f"Found {len(intersections)} intersections")
    
    # Visualize the vascular system and the intersections, and save CSV
    fig = vascular_system.visualize_intersection(plane_params, intersections, save_csv=True)
    plt.show()
    
    # Create a rotating animation of the 3D vascular system
    print("Creating animation of the 3D vascular system...")
    create_animation(vascular_system, plane_params)

def compare_different_seeds():
    """Compare vascular systems with different random seeds"""
    fig = plt.figure(figsize=(15, 10))
    
    for i, seed in enumerate([42, 123, 456, 789]):
        ax = fig.add_subplot(2, 2, i+1, projection='3d')
        
        # Create vascular system with this seed
        vascular_system = VascularSystem(
            start_point=[0, 0, 0],
            direction=[0, 0, 1],
            initial_length=10,
            initial_radius=2.0,
            min_branching_angle=20,
            max_branching_angle=45,
            radius_ratio=0.7,
            length_ratio=0.8,
            max_depth=9,
            min_radius=0.15,
            random_seed=seed,
            progressive_growth=False
        )
        
        # Visualize
        vascular_system.visualize_3d(ax)
        ax.set_title(f"Seed: {seed}")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Run the main function
    main()
    
    # Uncomment to compare different seeds (very slow)
    # compare_different_seeds()