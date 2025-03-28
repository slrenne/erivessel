import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
from dataclasses import dataclass
import math
from matplotlib.patches import Ellipse

@dataclass
class Vessel:
    start: np.ndarray  # 3D point where vessel starts
    end: np.ndarray    # 3D point where vessel ends
    radius: float      # vessel radius
    parent: object = None  # parent vessel

class VesselReconstructor:
    def __init__(self, slice_intersection, 
                 min_branching_angle=20, 
                 max_branching_angle=45, 
                 radius_ratio=0.7, 
                 length_ratio=0.8, 
                 max_depth=4, 
                 min_radius=0.1, 
                 random_seed=None):
        """
        Reconstruct 3D vascular system from a single 2D slice intersection
        
        Parameters:
        - slice_intersection: A single ellipse/circle intersection from the 2D slice
        - Other parameters similar to the original VascularSystem class
        """
        # Set random seed for reproducibility
        if random_seed is not None:
            random.seed(random_seed)
            np.random.seed(random_seed)
        
        self.slice_intersection = slice_intersection
        self.min_branching_angle = np.radians(min_branching_angle)
        self.max_branching_angle = np.radians(max_branching_angle)
        self.radius_ratio = radius_ratio
        self.length_ratio = length_ratio
        self.max_depth = max_depth
        self.min_radius = min_radius
        
        # Vessels will be stored here
        self.vessels = []
        self.vessels_by_generation = {}
        
        # Reconstruct the 3D vessel system
        self._reconstruct_vessels()
    
    def _normalize(self, vector):
        """Normalize a vector to unit length"""
        norm = np.linalg.norm(vector)
        return vector / norm if norm != 0 else vector
    
    def _get_perpendicular_vector(self, vector):
        """Find a perpendicular vector to the given vector"""
        if abs(vector[0]) < abs(vector[1]):
            perpendicular = np.array([1, 0, 0])
        else:
            perpendicular = np.array([0, 1, 0])
        
        return self._normalize(np.cross(vector, perpendicular))
    
    def _rotate_vector(self, vector, axis, angle):
        """Rotate a vector around an axis by the given angle"""
        axis = self._normalize(axis)
        k_cross_v = np.cross(axis, vector)
        return (vector * np.cos(angle) + 
                k_cross_v * np.sin(angle) + 
                axis * np.dot(axis, vector) * (1 - np.cos(angle)))
    
    def _reconstruct_vessels(self):
        """Reconstruct 3D vessels based on 2D slice intersection"""
        # Use the slice intersection center as the start point
        start_point = self.slice_intersection['center']
        
        # Use major axis as a basis for initial vessel length and radius
        initial_radius = self.slice_intersection['major_axis'] / 2
        initial_length = initial_radius * 3
        
        # Choose an initial 3D direction 
        # We'll use the plane normal and create a direction that's not parallel
        plane_normal = self.slice_intersection.get('normal', np.array([0, 0, 1]))
        
        # Create a base direction that's not parallel to the plane normal
        base_direction = np.array([1, 0, 0]) if not np.array_equal(plane_normal, [1, 0, 0]) else np.array([0, 1, 0])
        direction = self._normalize(base_direction - np.dot(base_direction, plane_normal) * plane_normal)
        
        # Create the root vessel
        end_point = start_point + direction * initial_length
        root_vessel = Vessel(
            start=start_point,
            end=end_point,
            radius=initial_radius
        )
        
        # Store the root vessel
        self.vessels.append(root_vessel)
        self.vessels_by_generation[0] = [root_vessel]
        
        # Generate branches from the root
        self._generate_vessels(root_vessel, 1)
    
    def _generate_vessels(self, parent_vessel, depth):
        """Recursively generate vessels with bifurcations"""
        if depth >= self.max_depth or parent_vessel.radius <= self.min_radius:
            return
        
        # Calculate child vessel parameters
        parent_direction = self._normalize(parent_vessel.end - parent_vessel.start)
        child_radius = parent_vessel.radius * self.radius_ratio
        child_length = np.linalg.norm(parent_vessel.end - parent_vessel.start) * self.length_ratio
        
        # Get a perpendicular axis for rotation
        rotation_axis = self._get_perpendicular_vector(parent_direction)
        
        # Increase branching factor as we go deeper
        branching_factor = min(4, 2 + int(depth/2))
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
    
    def visualize_3d(self, ax=None):
        """Visualize the reconstructed 3D vascular system"""
        if ax is None:
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(111, projection='3d')
        
        # Plot each vessel as a line segment
        for vessel in self.vessels:
            ax.plot(
                [vessel.start[0], vessel.end[0]],
                [vessel.start[1], vessel.end[1]],
                [vessel.start[2], vessel.end[2]],
                'r-', linewidth=vessel.radius * 5
            )
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Reconstructed 3D Vascular System')
        
        return ax

def generate_sample_slice_intersections(num_vessels=5):
    """
    Generate multiple 2D slice intersections representing vessel cross-sections
    
    Parameters:
    - num_vessels: Number of vessel intersections to generate
    
    Returns:
    - List of vessel intersection dictionaries
    """
    np.random.seed(42)
    
    intersections = []
    
    for _ in range(num_vessels):
        # Randomize center position
        center = np.array([
            np.random.uniform(-2, 2),  # x range
            np.random.uniform(-2, 2),  # y range
            0  # z-axis remains at 0 for 2D slice
        ])
        
        # Randomize vessel dimensions (with constraints)
        major_axis = np.random.uniform(1.0, 3.0)  # Diameter range
        minor_axis = major_axis * np.random.uniform(0.6, 1.0)  # Slight variation from circular
        
        # Randomize rotation angle
        rotation_angle = np.random.uniform(0, np.pi)
        
        # Randomize plane normal (with slight variation)
        normal = np.array([
            np.random.uniform(-0.2, 0.2),
            np.random.uniform(-0.2, 0.2),
            1  # Predominantly along z-axis
        ])
        normal = normal / np.linalg.norm(normal)
        
        intersection = {
            'center': center,
            'major_axis': major_axis,
            'minor_axis': minor_axis,
            'angle': rotation_angle,
            'normal': normal
        }
        
        intersections.append(intersection)
    
    return intersections

def main():
    # Generate multiple 2D slice intersections
    slice_intersections = generate_sample_slice_intersections(num_vessels=5)
    
    # Visualize original 2D slice intersections
    fig = plt.figure(figsize=(16, 6))
    
    # 2D slice plot
    ax1 = fig.add_subplot(121)
    ax1.set_title('Original 2D Slice Intersections')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_xlim(-3, 3)
    ax1.set_ylim(-3, 3)
    ax1.grid(True)
    ax1.set_aspect('equal')
    
    # Plot each intersection as an ellipse
    for intersection in slice_intersections:
        ellipse = Ellipse(
            xy=(intersection['center'][0], intersection['center'][1]), 
            width=intersection['major_axis'], 
            height=intersection['minor_axis'], 
            angle=np.degrees(intersection['angle']),
            fill=False, 
            color='red',
            alpha=0.7
        )
        ax1.add_artist(ellipse)
    
    # Reconstruct 3D vascular system from the slices
    # We'll use the first slice as the primary reconstruction point
    reconstructor = VesselReconstructor(
        slice_intersections[0],  # Primary slice
        min_branching_angle=20,
        max_branching_angle=45,
        radius_ratio=0.7,
        length_ratio=0.8,
        max_depth=4,
        random_seed=42
    )
    
    # 3D reconstruction plot
    ax2 = fig.add_subplot(122, projection='3d')
    reconstructor.visualize_3d(ax2)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()