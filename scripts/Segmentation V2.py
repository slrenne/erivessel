import cv2
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import ndimage
from scipy.spatial.distance import cdist


# ========================================
# CONFIGURATION PARAMETERS
# ========================================

class Config:
    """Configuration class for all segmentation parameters"""
    
    # Path settings
    INPUT_PATH = r"C:\Users\cirok\Desktop\Fix Segmentation\Input"
    OUTPUT_DIR = "Output_Segmentation"
    
    # File extensions to process
    SUPPORTED_EXTENSIONS = ['.jpeg', '.jpg', '.png', '.tiff', '.bmp']
    
    # Hole filling parameters
    HOLE_FILLING_KERNEL_SIZE = 15
    HOLE_FILLING_ITERATIONS = 2
    
    # Contour connection parameters
    MAX_CONNECTION_DISTANCE = 50
    MIN_ALIGNMENT_THRESHOLD = 0.3
    CONNECTION_LINE_THICKNESS = 3
    MIN_CONTOUR_AREA = 10
    
    # Progressive closing parameters
    MIN_KERNEL_SIZE = 3
    MAX_KERNEL_SIZE = 35
    KERNEL_STEP = 2
    MAX_AREA_INCREASE = 0.5  # 50% maximum area increase
    DILATION_ITERATIONS = 3
    
    # Component filtering parameters
    MIN_COMPONENT_AREA = 50
    CONTOUR_THICKNESS = 2
    CONTOUR_COLOR = (0, 255, 0)  # Green in BGR
    
    # Endpoint detection parameters
    ENDPOINT_SAMPLE_SIZE = 5


# ========================================
# UTILITY FUNCTIONS
# ========================================

def setup_output_directory(output_dir):
    """Create output directory if it doesn't exist"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")


def get_image_files(path, extensions):
    """Get all image files with specified extensions from the input directory"""
    image_files = []
    for ext in extensions:
        image_files.extend([f for f in os.listdir(path) if f.lower().endswith(ext.lower())])
    return image_files


def create_image_output_directory(output_dir, filename):
    """Create subdirectory for individual image results"""
    image_output_dir = os.path.join(output_dir, f"{filename}_SEGMENTATION")
    if not os.path.exists(image_output_dir):
        os.makedirs(image_output_dir)
    return image_output_dir


# ========================================
# HOLE FILLING FUNCTIONS
# ========================================

def fill_holes_comprehensive(binary, kernel_size=None, iterations=None):
    """
    Comprehensive hole filling using multiple approaches
    
    Args:
        binary: Binary image to process
        kernel_size: Size of morphological kernel (default from config)
        iterations: Number of iterations (default from config)
    """
    if kernel_size is None:
        kernel_size = Config.HOLE_FILLING_KERNEL_SIZE
    if iterations is None:
        iterations = Config.HOLE_FILLING_ITERATIONS
    
    result = binary.copy()
    h, w = binary.shape
    
    # Method 1: Flood fill from borders
    filled = binary.copy()
    mask = np.zeros((h + 2, w + 2), dtype=np.uint8)
    cv2.floodFill(filled, mask, (0, 0), 255)
    filled_inv = cv2.bitwise_not(filled)
    result = cv2.bitwise_or(result, filled_inv)
    
    # Method 2: Fill contours
    contours, _ = cv2.findContours(result, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    filled_contours = np.zeros_like(result)
    for contour in contours:
        cv2.fillPoly(filled_contours, [contour], 255)
    
    # Method 3: Morphological closing
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size))
    closed_holes = cv2.morphologyEx(result, cv2.MORPH_CLOSE, kernel, iterations=iterations)
    
    # Combine all methods
    final_result = cv2.bitwise_or(result, filled_contours)
    final_result = cv2.bitwise_or(final_result, closed_holes)
    
    return final_result


def advanced_hole_filling(binary):
    """Advanced hole filling using contour hierarchy and ndimage"""
    result = binary.copy()
    
    # Fill holes using contour hierarchy
    contours, hierarchy = cv2.findContours(result, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    filled = np.zeros_like(result)
    
    # Fill all external contours
    for i, contour in enumerate(contours):
        if hierarchy[0][i][3] == -1:  # External contour
            cv2.drawContours(filled, [contour], -1, 255, -1)
    
    # Additional hole filling using ndimage
    filled_ndimage = ndimage.binary_fill_holes(filled).astype(np.uint8) * 255
    
    # Combine with original
    result = cv2.bitwise_or(result, filled_ndimage)
    
    return result


# ========================================
# CONTOUR CONNECTION FUNCTIONS
# ========================================

def get_endpoints_and_orientation(contour, num_points=None):
    """
    Get endpoints and local orientation of a contour
    
    Args:
        contour: Input contour
        num_points: Number of points to sample (default from config)
    """
    if num_points is None:
        num_points = Config.ENDPOINT_SAMPLE_SIZE
    
    if len(contour) < num_points:
        return None, None, None, None
    
    # Get start and end points
    start_points = contour[:num_points].reshape(-1, 2)
    end_points = contour[-num_points:].reshape(-1, 2)
    
    start_center = np.mean(start_points, axis=0)
    end_center = np.mean(end_points, axis=0)
    
    # Calculate orientations
    start_direction = start_points[-1] - start_points[0]
    end_direction = end_points[-1] - end_points[0]
    
    # Normalize
    start_direction = start_direction / (np.linalg.norm(start_direction) + 1e-6)
    end_direction = end_direction / (np.linalg.norm(end_direction) + 1e-6)
    
    return start_center, end_center, start_direction, end_direction


def should_connect_contours(c1, c2, max_distance=None, min_alignment=None):
    """
    Determine if two contours should be connected
    
    Args:
        c1, c2: Contours to evaluate
        max_distance: Maximum connection distance (default from config)
        min_alignment: Minimum alignment threshold (default from config)
    """
    if max_distance is None:
        max_distance = Config.MAX_CONNECTION_DISTANCE
    if min_alignment is None:
        min_alignment = Config.MIN_ALIGNMENT_THRESHOLD
    
    # Get endpoint information
    start1, end1, start_dir1, end_dir1 = get_endpoints_and_orientation(c1)
    start2, end2, start_dir2, end_dir2 = get_endpoints_and_orientation(c2)
    
    if any(x is None for x in [start1, end1, start_dir1, end_dir1, start2, end2, start_dir2, end_dir2]):
        return False
    
    # Check all possible connections
    connections = [
        (end1, start2, -end_dir1, start_dir2),
        (end1, end2, -end_dir1, -end_dir2),
        (start1, start2, start_dir1, start_dir2),
        (start1, end2, start_dir1, -end_dir2)
    ]
    
    for p1, p2, dir1, dir2 in connections:
        distance = np.linalg.norm(p1 - p2)
        if distance <= max_distance:
            alignment = np.dot(dir1, dir2)
            if alignment > min_alignment:
                return True
    
    return False


def connect_nearby_contours(binary, max_distance=None, min_alignment=None):
    """
    Connect nearby contours that are likely part of the same vessel
    
    Args:
        binary: Binary image
        max_distance: Maximum connection distance (default from config)
        min_alignment: Minimum alignment threshold (default from config)
    """
    if max_distance is None:
        max_distance = Config.MAX_CONNECTION_DISTANCE
    if min_alignment is None:
        min_alignment = Config.MIN_ALIGNMENT_THRESHOLD
    
    # Find and filter contours
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    contours = [c for c in contours if cv2.contourArea(c) > Config.MIN_CONTOUR_AREA]
    
    result = binary.copy()
    
    # Check each pair of contours
    for i in range(len(contours)):
        for j in range(i + 1, len(contours)):
            if should_connect_contours(contours[i], contours[j], max_distance, min_alignment):
                # Connect contours with a line
                points1 = contours[i].reshape(-1, 2)
                points2 = contours[j].reshape(-1, 2)
                
                # Find closest points
                distances = cdist(points1, points2)
                min_idx = np.unravel_index(np.argmin(distances), distances.shape)
                
                pt1 = tuple(points1[min_idx[0]])
                pt2 = tuple(points2[min_idx[1]])
                
                # Draw connection line
                cv2.line(result, pt1, pt2, 255, Config.CONNECTION_LINE_THICKNESS)
    
    return result


# ========================================
# PROGRESSIVE CLOSING FUNCTIONS
# ========================================

def enhanced_progressive_closing(binary, max_kernel_size=None):
    """
    Enhanced progressive closing with precision control
    
    Args:
        binary: Binary image to process
        max_kernel_size: Maximum kernel size (default from config)
    """
    if max_kernel_size is None:
        max_kernel_size = Config.MAX_KERNEL_SIZE
    
    result = binary.copy()
    
    # Progressive closing with increasing kernel sizes
    for kernel_size in range(Config.MIN_KERNEL_SIZE, max_kernel_size + 1, Config.KERNEL_STEP):
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size))
        
        # More conservative iterations based on kernel size
        iterations = 1  # Always use 1 iteration for precision
        closed = cv2.morphologyEx(result, cv2.MORPH_CLOSE, kernel, iterations=iterations)
        
        # Stricter validation - check both area and shape changes
        original_area = cv2.countNonZero(result)
        new_area = cv2.countNonZero(closed)
        area_increase = (new_area - original_area) / original_area if original_area > 0 else 0
        
        # More stringent area control for larger kernels
        max_allowed_increase = Config.MAX_AREA_INCREASE
        if kernel_size > 10:
            max_allowed_increase = Config.MAX_AREA_INCREASE * 0.5  # Even stricter for large kernels
        
        if area_increase < max_allowed_increase:
            # Additional validation: check if we're merging separate components
            original_components = cv2.connectedComponents(result)[0]
            new_components = cv2.connectedComponents(closed)[0]
            
            # If we're significantly reducing the number of components, be more careful
            if new_components >= original_components * 0.8:  # Allow max 20% reduction in components
                # Create boundary validation with reduced dilation
                original_dilated = cv2.dilate(binary, kernel, iterations=max(1, Config.DILATION_ITERATIONS // 2))
                valid_changes = cv2.bitwise_and(closed, original_dilated)
                result = cv2.bitwise_or(result, valid_changes)
    
    # Apply single round of conservative hole filling
    result = fill_holes_comprehensive(result)
    
    return result


# ========================================
# MAIN SEGMENTATION FUNCTION
# ========================================

def ultimate_vessel_segmentation(binary):
    """
    Ultimate vessel segmentation with comprehensive processing pipeline
    
    Args:
        binary: Binary image to process
        
    Returns:
        Processed binary image with filled vessels
    """
    # Step 1: Connect nearby contours
    connected = connect_nearby_contours(binary)
    
    # Step 2: Progressive closing
    closed = enhanced_progressive_closing(connected)
    
    # Step 3: Advanced hole filling
    filled = advanced_hole_filling(closed)
    
    # Step 4: Remove small components
    num_labels, labels = cv2.connectedComponents(filled)
    final_result = np.zeros_like(filled)
    
    for label in range(1, num_labels):
        component_mask = (labels == label).astype(np.uint8) * 255
        area = cv2.countNonZero(component_mask)
        
        if area > Config.MIN_COMPONENT_AREA:
            final_result = cv2.bitwise_or(final_result, component_mask)
    
    # Step 5: Final hole filling
    final_result = fill_holes_comprehensive(final_result)
    
    return final_result


# ========================================
# IMAGE PROCESSING PIPELINE
# ========================================

def process_image(image_path, output_dir):
    """
    Process a single image through the vessel segmentation pipeline
    
    Args:
        image_path: Path to the input image
        output_dir: Directory to save results
    """
    # Load image
    original = cv2.imread(image_path)
    if original is None:
        print(f"Error: Could not read image {image_path}")
        return False
    
    # Get filename without extension
    filename = os.path.splitext(os.path.basename(image_path))[0]
    print(f"Processing image: {filename}")
    
    # Create output directory for this image
    image_output_dir = create_image_output_directory(output_dir, filename)
    
    # Convert to HSV and extract H channel
    hsv_image = cv2.cvtColor(original, cv2.COLOR_BGR2HSV)
    h_channel = cv2.split(hsv_image)[0]
    h_channel = 255 - h_channel
    
    # Apply OTSU thresholding
    _, binary = cv2.threshold(h_channel, 0, 255, cv2.THRESH_OTSU)
    binary = binary.astype(np.uint8)
    
    # Apply enhanced progressive closing
    result = enhanced_progressive_closing(binary)
    
    # Save binary result
    binary_output_path = os.path.join(image_output_dir, f"{filename}_Enhanced_Progressive_Filled.png")
    plt.imsave(binary_output_path, result, cmap="gray")
    
    # Create and save overlay
    contours, _ = cv2.findContours(result, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    overlay_image = original.copy()
    
    for contour in contours:
        cv2.drawContours(overlay_image, [contour], -1, Config.CONTOUR_COLOR, Config.CONTOUR_THICKNESS)
    
    overlay_output_path = os.path.join(image_output_dir, f"{filename}_Vessels_Enhanced_Progressive_Filled.jpg")
    cv2.imwrite(overlay_output_path, overlay_image)
    
    return True


# ========================================
# MAIN EXECUTION
# ========================================

def main():
    """Main execution function"""
    print("=" * 60)
    print("VESSEL SEGMENTATION PIPELINE")
    print("=" * 60)
    
    # Setup output directory
    setup_output_directory(Config.OUTPUT_DIR)
    
    # Get image files
    image_files = get_image_files(Config.INPUT_PATH, Config.SUPPORTED_EXTENSIONS)
    
    if not image_files:
        print(f"No supported image files found in {Config.INPUT_PATH}")
        return
    
    print(f"Found {len(image_files)} image(s) to process")
    
    # Process each image
    processed_count = 0
    for image_file in image_files:
        image_path = os.path.join(Config.INPUT_PATH, image_file)
        if process_image(image_path, Config.OUTPUT_DIR):
            processed_count += 1
    
    # Summary
    print("=" * 60)
    print(f"Processing complete!")
    print(f"Successfully processed: {processed_count}/{len(image_files)} images")
    print(f"Results saved in: {Config.OUTPUT_DIR}")
    print("Applied Enhanced Progressive Closing with comprehensive hole filling")
    print("=" * 60)


if __name__ == "__main__":
    main()