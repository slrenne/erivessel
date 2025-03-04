import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse
import random
from dataclasses import dataclass
import math
from matplotlib.animation import FuncAnimation

@dataclass
class Vessel:
    start: np.ndarray  # 3D point where vessel starts
    end: np.ndarray    # 3D point where vessel ends
    radius: float      # vessel radius
    parent: object = None  # parent vessel

class VesselReconstructor:
    def __init__(self, plane_params=None, random_seed=None):
        """
        Initialize a vessel reconstructor that builds 3D vessels from 2D slice data
        
        Parameters:
        - plane_params: Optional tuple (a, b, c, d) defining the cutting plane ax + by + cz + d = 0
        - random_seed: Optional seed for random number generation
        """
        # Set random seed if provided
        if random_seed is not None:
            random.seed(random_seed)
            np.random.seed(random_seed)
        
        # Set default plane if none provided
        if plane_params is None:
            # Default: horizontal plane at z=0
            self.plane_params = (0, 0, 1, 0)  # z = 0 plane
        else:
            self.plane_params = plane_params
            
        # Normalize the normal vector of the plane
        a, b, c, d = self.plane_params
        normal = np.array([a, b, c])
        self.normal = self._normalize(normal)
        
        # Create basis vectors for the plane
        self.u = self._get_perpendicular_vector(self.normal)
        self.v = np.cross(self.normal, self.u)
        
        # Find a point on the plane (origin for the 2D coordinate system)
        if c != 0:
            self.origin = np.array([0, 0, -d/c])
        elif b != 0:
            self.origin = np.array([0, -d/b, 0])
        else:
            self.origin = np.array([-d/a, 0, 0])
            
        # List to store 3D reconstructed vessels
        self.vessels = []
        
        # List to store 2D ellipses
        self.ellipses = []
    
    def _normalize(self, vector):
        """Normalize a vector to unit length"""
        norm = np.linalg.norm(vector)
        if norm == 0:
            return vector
        return vector / norm
    
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
    
    def generate_random_slice(self, num_vessels=10, area_size=20, min_radius=0.2, max_radius=2.0,
                              min_eccentricity=0.0, max_eccentricity=0.8):
        """
        Generate a random 2D slice with elliptical vessel cross-sections
        
        Parameters:
        - num_vessels: Number of vessel cross-sections to generate
        - area_size: Size of the 2D area for generation
        - min_radius: Minimum vessel radius
        - max_radius: Maximum vessel radius
        - min_eccentricity: Minimum eccentricity (0=circle, approaching 1=very elongated)
        - max_eccentricity: Maximum eccentricity
        
        Returns:
        - List of dictionaries containing ellipse parameters
        """
        self.ellipses = []
        
        for _ in range(num_vessels):
            # Generate random position
            x = random.uniform(-area_size/2, area_size/2)
            y = random.uniform(-area_size/2, area_size/2)
            
            # Generate random radius and eccentricity
            minor_axis = random.uniform(min_radius, max_radius)
            eccentricity = random.uniform(min_eccentricity, max_eccentricity)
            major_axis = minor_axis / (1 - eccentricity) if eccentricity < 1 else minor_axis * 5
            
            # Random rotation angle
            angle = random.uniform(0, 2 * np.pi)
            
            # Store the ellipse data
            ellipse = {
                'center': np.array([x, y]),
                'major_axis': major_axis,
                'minor_axis': minor_axis,
                'angle': angle
            }
            self.ellipses.append(ellipse)
        
        return self.ellipses
    
    def _2d_to_3d_point(self, point_2d):
        """Convert a 2D point in the plane to a 3D point"""
        return self.origin + point_2d[0] * self.u + point_2d[1] * self.v
    
    def reconstruct_3d_vessels(self, vessel_length=10.0, min_angle=20, max_angle=70, 
                               bidirectional=True):
        """
        Reconstruct 3D vessels from 2D elliptical cross-sections
        
        Parameters:
        - vessel_length: Length of the reconstructed vessels
        - min_angle: Minimum angle (in degrees) between vessel and normal
        - max_angle: Maximum angle (in degrees) between vessel and normal
        - bidirectional: If True, create vessels in both directions from the slice
        
        Returns:
        - List of Vessel objects
        """
        self.vessels = []
        
        for ellipse in self.ellipses:
            # Convert 2D center to 3D point
            center_3d = self._2d_to_3d_point(ellipse['center'])
            
            # Get the ellipse parameters
            major_axis = ellipse['major_axis']
            minor_axis = ellipse['minor_axis']
            angle = ellipse['angle']
            
            # The minor axis of the ellipse corresponds to the vessel radius
            vessel_radius = minor_axis
            
            # Calculate the vessel direction based on the ellipse properties
            # For a circle (major_axis = minor_axis), the vessel is perpendicular to the plane
            # For an ellipse, the vessel's angle depends on the ratio of major_axis to minor_axis
            
            if np.isclose(major_axis, minor_axis):
                # For a circle, the vessel is perpendicular to the plane
                vessel_direction = self.normal.copy()
            else:
                # For an ellipse, calculate the vessel angle from the major/minor axis ratio
                # sin_angle = minor_axis / major_axis (from the original code's transformation)
                sin_angle = minor_axis / major_axis
                cos_angle = np.sqrt(1 - sin_angle**2)
                
                # The vessel direction is at an angle from the normal
                # First, get the in-plane direction aligned with the ellipse's major axis
                in_plane_dir = self._rotate_vector(self.u, self.normal, angle)
                
                # Then, rotate away from normal by the calculated angle
                vessel_direction = cos_angle * self.normal + sin_angle * in_plane_dir
                vessel_direction = self._normalize(vessel_direction)
                
                # Ensure the angle is within the specified range
                current_angle_deg = np.degrees(np.arccos(np.dot(vessel_direction, self.normal)))
                if not (min_angle <= current_angle_deg <= max_angle):
                    # Calculate a new angle within range
                    new_angle_rad = np.radians(random.uniform(min_angle, max_angle))
                    sin_angle = np.sin(new_angle_rad)
                    cos_angle = np.cos(new_angle_rad)
                    
                    # Recalculate vessel direction
                    vessel_direction = cos_angle * self.normal + sin_angle * in_plane_dir
                    vessel_direction = self._normalize(vessel_direction)
            
            # Calculate vessel endpoints
            vessel_half_length = vessel_length / 2 if bidirectional else vessel_length
            
            if bidirectional:
                start_3d = center_3d - vessel_direction * vessel_half_length
                end_3d = center_3d + vessel_direction * vessel_half_length
            else:
                start_3d = center_3d
                end_3d = center_3d + vessel_direction * vessel_length
            
            # Create the vessel
            vessel = Vessel(
                start=start_3d,
                end=end_3d,
                radius=vessel_radius
            )
            
            self.vessels.append(vessel)
        
        return self.vessels
    
    def visualize_2d_slice(self, ax=None):
        """Visualize the 2D slice with vessel cross-sections"""
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 10))
            standalone = True
        else:
            standalone = False
        
        # Plot each ellipse
        for ellipse_data in self.ellipses:
            center = ellipse_data['center']
            major_axis = ellipse_data['major_axis']
            minor_axis = ellipse_data['minor_axis']
            angle_rad = ellipse_data['angle']
            angle_deg = np.degrees(angle_rad)
            
            ellipse = Ellipse(
                xy=center,
                width=2 * major_axis,
                height=2 * minor_axis,
                angle=angle_deg,
                edgecolor='red',
                facecolor='none',
                linewidth=2
            )
            
            ax.add_patch(ellipse)
        
        # Set aspect ratio and labels
        ax.set_aspect('equal')
        
        # Calculate plot limits
        if self.ellipses:
            centers = np.array([e['center'] for e in self.ellipses])
            max_radius = max(max(e['major_axis'], e['minor_axis']) for e in self.ellipses)
            
            min_x, max_x = np.min(centers[:, 0]), np.max(centers[:, 0])
            min_y, max_y = np.min(centers[:, 1]), np.max(centers[:, 1])
            
            margin = max_radius * 2
            ax.set_xlim(min_x - margin, max_x + margin)
            ax.set_ylim(min_y - margin, max_y + margin)
        else:
            ax.set_xlim(-10, 10)
            ax.set_ylim(-10, 10)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('2D Vessel Cross-Sections')
        ax.grid(True)
        
        if standalone:
            plt.tight_layout()
            plt.show()
            return fig
        
        return ax
    
    def visualize_3d_vessels(self, ax=None, show_plane=True):
        """Visualize the reconstructed 3D vessels"""
        if ax is None:
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(111, projection='3d')
            standalone = True
        else:
            standalone = False
        
        # Plot each vessel as a line segment
        for vessel in self.vessels:
            ax.plot(
                [vessel.start[0], vessel.end[0]],
                [vessel.start[1], vessel.end[1]],
                [vessel.start[2], vessel.end[2]],
                'r-', linewidth=vessel.radius * 3  # Scale the linewidth
            )
        
        # Show the intersection plane if requested
        if show_plane:
            self._plot_plane(ax)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Reconstructed 3D Vessels')
        
        if standalone:
            plt.tight_layout()
            plt.show()
            return fig
        
        return ax
    
    def _plot_plane(self, ax):
        """Plot the intersection plane in 3D"""
        a, b, c, d = self.plane_params
        
        # Find the limits of our vascular system
        if self.vessels:
            xs = [v.start[0] for v in self.vessels] + [v.end[0] for v in self.vessels]
            ys = [v.start[1] for v in self.vessels] + [v.end[1] for v in self.vessels]
            zs = [v.start[2] for v in self.vessels] + [v.end[2] for v in self.vessels]
            
            min_x, max_x = min(xs), max(xs)
            min_y, max_y = min(ys), max(ys)
            margin = 0.2 * max(max_x - min_x, max_y - min_y)
            
            if margin < 1:  # Ensure minimum margin
                margin = 5
        else:
            # Default if no vessels
            min_x, max_x = -10, 10
            min_y, max_y = -10, 10
            margin = 5
        
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
            if self.vessels:
                min_z, max_z = min(zs), max(zs)
            else:
                min_z, max_z = -10, 10
                
            z_range = np.linspace(min_z - margin, max_z + margin, 10)
            
            if a != 0:  # Plane is of form ax + by + d = 0 (c=0)
                yy, zz = np.meshgrid(y_range, z_range)
                xx = (-b * yy - d) / a
                ax.plot_surface(xx, yy, zz, alpha=0.3, color='blue')
            elif b != 0:  # Plane is of form by + d = 0 (a=0, c=0)
                xx, zz = np.meshgrid(x_range, z_range)
                yy = np.ones_like(xx) * (-d / b)
                ax.plot_surface(xx, yy, zz, alpha=0.3, color='blue')
    
    def visualize_combined(self):
        """Show both 2D slice and 3D reconstruction side by side"""
        fig = plt.figure(figsize=(20, 10))
        
        # 2D slice visualization
        ax1 = fig.add_subplot(121)
        self.visualize_2d_slice(ax1)
        
        # 3D vessel visualization
        ax2 = fig.add_subplot(122, projection='3d')
        self.visualize_3d_vessels(ax2)
        
        plt.tight_layout()
        plt.show()
        return fig
    
    def create_animation(self, save_animation=False):
        """Create an animation of rotating the 3D vascular system"""
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Visualize the system
        self.visualize_3d_vessels(ax)
        
        # Set up the animation
        def update(frame):
            ax.view_init(elev=20, azim=frame)
            return fig,
        
        # Create the animation
        ani = FuncAnimation(fig, update, frames=range(0, 360, 2), interval=50, blit=True)
        
        # Save if requested
        if save_animation:
            ani.save('reconstructed_vascular_system_rotation.gif', writer='pillow', fps=20)
        
        plt.tight_layout()
        plt.show()
        
        return ani

def random_plane_params(min_val=-1, max_val=1):
    """Generate random plane parameters"""
    a = random.uniform(min_val, max_val)
    b = random.uniform(min_val, max_val)
    c = random.uniform(min_val, max_val)
    
    # Normalize the normal vector
    norm = np.sqrt(a**2 + b**2 + c**2)
    a, b, c = a/norm, b/norm, c/norm
    
    # Random offset (d parameter)
    d = random.uniform(min_val, max_val)
    
    return (a, b, c, d)

def main():
    # Set random seed for reproducibility
    seed = 42
    random.seed(seed)
    np.random.seed(seed)
    
    # Create a vessel reconstructor (optionally with a random plane)
    plane_params = random_plane_params()
    print(f"Plane equation: {plane_params[0]:.3f}x + {plane_params[1]:.3f}y + {plane_params[2]:.3f}z + {plane_params[3]:.3f} = 0")
    
    reconstructor = VesselReconstructor(plane_params=plane_params, random_seed=seed)
    
    # Generate a random 2D slice with vessel cross-sections
    print("Generating random 2D vessel cross-sections...")
    reconstructor.generate_random_slice(
        num_vessels=15,
        area_size=20,
        min_radius=0.3,
        max_radius=1.5,
        min_eccentricity=0.0,  # 0 = circle
        max_eccentricity=0.8   # Close to 1 = very elongated ellipse
    )
    
    # Reconstruct 3D vessels from the 2D slice
    print("Reconstructing 3D vessels from the 2D slice...")
    reconstructor.reconstruct_3d_vessels(
        vessel_length=15.0,
        min_angle=20,   # Minimum angle between vessel and normal
        max_angle=70,   # Maximum angle between vessel and normal
        bidirectional=True  # Vessels extend in both directions from the slice
    )
    
    # Visualize both 2D slice and 3D reconstruction
    print("Visualizing the results...")
    reconstructor.visualize_combined()
    
    # Create a rotating animation of the 3D reconstruction
    print("Creating animation of the 3D reconstruction...")
    reconstructor.create_animation()

if __name__ == "__main__":
    main()