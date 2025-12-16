import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian
from skimage.segmentation import active_contour
from skimage.measure import label, regionprops
from skimage import measure
from scipy.ndimage import binary_fill_holes
import pandas as pd
import os

def process_single_image(image_path, output_dir=None, show_plots=True, save_plots=True,
                         # Morphological parameters for vessel closure
                         blur_kernel=5, 
                         closing_kernel=5, 
                         dilation_kernel=21, dilation_iterations=1,
                         erosion_kernel=3, erosion_iterations=2,
                         cleanup_kernel=5):
    """
    Process a single image using CMYK segmentation and extract measurements with filled lumens
    
    Parameters:
    image_path: path to the input image
    output_dir: directory to save results (if None, creates one based on image name)
    show_plots: whether to display plots
    save_plots: whether to save plots and results
    
    Morphological parameters for vessel closure:
    blur_kernel: size of Gaussian blur kernel (higher = more smoothing)
    closing_kernel: size of morphological closing kernel (higher = closes larger gaps)
    dilation_kernel: size of dilation kernel (higher = connects more distant structures)
    dilation_iterations: number of dilation iterations (higher = more aggressive connection)
    erosion_kernel: size of erosion kernel (higher = more refinement)
    erosion_iterations: number of erosion iterations (higher = more refinement)
    cleanup_kernel: size of final cleanup kernel (higher = smoother result)
    """
    
    # Create output directory if saving is enabled
    if save_plots:
        if output_dir is None:
            base_name = os.path.splitext(os.path.basename(image_path))[0]
            output_dir = f"{base_name}_CMYK_ANALYSIS"
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    print(f"Processing image: {image_path}")
    
    # Step 1: Load and convert image
    img_bgr = cv2.imread(image_path)
    if img_bgr is None:
        print(f"Error: Could not load image {image_path}")
        return None
    
    # Make float and divide by 255 to give BGRdash
    bgrdash = img_bgr.astype(np.float64)/255.

    # Calculate K as (1 - whatever is biggest out of Rdash, Gdash, Bdash)
    K = 1 - np.max(bgrdash, axis=2)

    # Calculate C
    C = (1-bgrdash[...,2] - K)/(1-K)

    # Calculate M
    M = (1-bgrdash[...,1] - K)/(1-K)

    # Calculate Y
    Y = (1-bgrdash[...,0] - K)/(1-K)

    # Combine 4 channels into single image and re-scale back up to uint8
    CMYK = (np.dstack((C,M,Y,K)) * 255).astype(np.uint8)
    
    # Apply Gaussian blur to help connect nearby structures
    Y_channel = CMYK[:,:,2]
    Y_blurred = cv2.GaussianBlur(Y_channel, (blur_kernel, blur_kernel), 0)
    
    # Apply Otsu thresholding on blurred Yellow channel
    _, binary = cv2.threshold(Y_blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    # Morphological operations for better vessel closure
    # First, apply morphological closing to connect nearby vessel segments
    closing_kernel_elem = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (closing_kernel, closing_kernel))
    closed = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, closing_kernel_elem)
    
    # Then dilate to further connect structures
    dilation_kernel_elem = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (dilation_kernel, dilation_kernel))
    dilated = cv2.dilate(closed, dilation_kernel_elem, iterations=dilation_iterations)
    
    # Finally erode to refine
    erosion_kernel_elem = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (erosion_kernel, erosion_kernel))
    eroded = cv2.erode(dilated, erosion_kernel_elem, iterations=erosion_iterations)

    # NEW: Fill holes to create solid lumens
    filled = binary_fill_holes(eroded).astype(np.uint8) * 255
    
    # Optional: Additional morphological operations to clean up the filled result
    # You can adjust these parameters based on your needs
    kernel_cleanup = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (cleanup_kernel, cleanup_kernel))
    filled_cleaned = cv2.morphologyEx(filled, cv2.MORPH_CLOSE, kernel_cleanup)
    filled_cleaned = cv2.morphologyEx(filled_cleaned, cv2.MORPH_OPEN, kernel_cleanup)

    # Find Contours on the filled image
    contours, _ = cv2.findContours(filled_cleaned, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Approximate each contour with tighter precision
    refined_contours = []
    for cnt in contours:
        epsilon = 0.0000001 * cv2.arcLength(cnt, True)  # Lower = tighter fit
        approx = cv2.approxPolyDP(cnt, epsilon, True)
        refined_contours.append(approx)
    
    # Create a copy of the original image for visualization
    contour_img = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB).copy()

    # Draw filled contours instead of just outlines
    cv2.drawContours(contour_img, refined_contours, -1, (180, 255, 0), -1)  # -1 thickness fills the contour
    
    # Also create an overlay version for better visualization
    overlay_img = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB).copy()
    mask = np.zeros(overlay_img.shape[:2], dtype=np.uint8)
    cv2.drawContours(mask, refined_contours, -1, 255, -1)
    
    # Create colored overlay
    overlay_img[mask > 0] = overlay_img[mask > 0] * 0.7 + np.array([180, 255, 0]) * 0.3

    # MEASUREMENT PART (updated to use filled image)
    # Define the properties to measure
    all_props = ["area", "axis_major_length", "axis_minor_length", "eccentricity", "orientation"]
    
    # Create labeled image for measurements using filled image
    label_img = label(filled_cleaned)
    
    # Calculate properties
    props = measure.regionprops_table(label_img, img_bgr, properties=all_props)
    
    # Apply area scaling (from original segmentation code)
    props["area"] = props["area"] * 0.25 ** 2
    
    # Create DataFrame
    data = pd.DataFrame(props)
    
    # Add additional statistics
    if len(data) > 0:
        print(f"Found {len(data)} objects")
        print(f"Total area: {data['area'].sum():.2f}")
        print(f"Mean area: {data['area'].mean():.2f}")
        print(f"Mean eccentricity: {data['eccentricity'].mean():.3f}")
    else:
        print("No objects detected")
    
    # Save results if requested
    if save_plots and output_dir:
        # Save measurements to CSV
        data.to_csv(f'{output_dir}/measurements.csv', index=True)
        
        # Save segmented image (now filled)
        plt.imsave(f"{output_dir}/segmented_filled.png", filled_cleaned, cmap="gray")
        
        # Save contoured image
        plt.imsave(f"{output_dir}/contoured_filled.png", contour_img)
        
        # Save overlay image
        plt.imsave(f"{output_dir}/overlay_filled.png", overlay_img)
        
        # Save the grid plot with additional steps
        fig, ax = plt.subplots(2, 4, figsize=(20, 10))

        ax[0,0].imshow(cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB))
        ax[0,0].set_title("Original Image")
        ax[0,0].axis("off")

        ax[0,1].imshow(Y_channel, cmap="gray")
        ax[0,1].set_title("Yellow Channel (Y)")
        ax[0,1].axis("off")

        ax[0,2].imshow(binary, cmap="gray")
        ax[0,2].set_title("Binarized (Otsu)")
        ax[0,2].axis("off")

        ax[0,3].imshow(dilated, cmap="gray")
        ax[0,3].set_title("Dilated (21x21 Kernel)")
        ax[0,3].axis("off")

        ax[1,0].imshow(eroded, cmap="gray")
        ax[1,0].set_title("Eroded (3x3 Kernel)")
        ax[1,0].axis("off")

        ax[1,1].imshow(filled_cleaned, cmap="gray")
        ax[1,1].set_title("Filled Lumens")
        ax[1,1].axis("off")

        ax[1,2].imshow(contour_img)
        ax[1,2].set_title("Filled Contours")
        ax[1,2].axis("off")

        ax[1,3].imshow(overlay_img)
        ax[1,3].set_title("Overlay")
        ax[1,3].axis("off")

        plt.tight_layout()
        plt.savefig(f"{output_dir}/processing_steps_filled.png", dpi=300)
        if show_plots:
            plt.show()
        plt.close()
        
        print(f"Results saved to: {output_dir}")
    
    elif show_plots:
        # Display plots
        fig, ax = plt.subplots(2, 4, figsize=(20, 10))

        ax[0,0].imshow(cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB))
        ax[0,0].set_title("Original Image")
        ax[0,0].axis("off")

        ax[0,1].imshow(Y_channel, cmap="gray")
        ax[0,1].set_title("Yellow Channel (Y)")
        ax[0,1].axis("off")

        ax[0,2].imshow(binary, cmap="gray")
        ax[0,2].set_title("Binarized (Otsu)")
        ax[0,2].axis("off")

        ax[0,3].imshow(dilated, cmap="gray")
        ax[0,3].set_title("Dilated (21x21 Kernel)")
        ax[0,3].axis("off")

        ax[1,0].imshow(eroded, cmap="gray")
        ax[1,0].set_title("Eroded (3x3 Kernel)")
        ax[1,0].axis("off")

        ax[1,1].imshow(filled_cleaned, cmap="gray")
        ax[1,1].set_title("Filled Lumens")
        ax[1,1].axis("off")

        ax[1,2].imshow(contour_img)
        ax[1,2].set_title("Filled Contours")
        ax[1,2].axis("off")

        ax[1,3].imshow(overlay_img)
        ax[1,3].set_title("Overlay")
        ax[1,3].axis("off")

        plt.tight_layout()
        plt.show()
        
        plt.figure(figsize=(10, 8))
        plt.imshow(overlay_img)
        plt.title("Final Filled Lumens Overlay")
        plt.axis("off")
        plt.show()
    
    return data

def process_folder(folder_path, output_base_dir=None, show_plots=False, save_plots=True):
    """
    Process all images in a folder using CMYK segmentation and extract measurements with filled lumens
    
    Parameters:
    folder_path: path to the folder containing images
    output_base_dir: base directory to save results (if None, creates one based on folder name)
    show_plots: whether to display plots for each image
    save_plots: whether to save plots and results
    """
    
    # Supported image extensions
    supported_extensions = ['.jpg', '.jpeg', '.png', '.bmp', '.tiff', '.tif']
    
    # Get all image files in the folder
    image_files = []
    for file in os.listdir(folder_path):
        if any(file.lower().endswith(ext) for ext in supported_extensions):
            image_files.append(os.path.join(folder_path, file))
    
    if not image_files:
        print(f"No supported image files found in {folder_path}")
        return
    
    print(f"Found {len(image_files)} image(s) to process")
    
    # Create base output directory if saving is enabled
    if save_plots:
        if output_base_dir is None:
            folder_name = os.path.basename(folder_path.rstrip('/\\'))
            output_base_dir = f"{folder_name}_CMYK_BATCH_ANALYSIS_FILLED"
        
        if not os.path.exists(output_base_dir):
            os.makedirs(output_base_dir)
    
    # Process each image
    all_measurements = []
    
    for i, image_path in enumerate(image_files):
        print(f"\n--- Processing image {i+1}/{len(image_files)} ---")
        
        # Create individual output directory for this image
        if save_plots:
            image_name = os.path.splitext(os.path.basename(image_path))[0]
            image_output_dir = os.path.join(output_base_dir, image_name)
        else:
            image_output_dir = None
        
        # Process the image
        measurements = process_single_image(image_path, image_output_dir, show_plots, save_plots)
        
        if measurements is not None and len(measurements) > 0:
            # Add image identifier to measurements
            measurements['image_name'] = os.path.basename(image_path)
            measurements['image_number'] = i + 1
            all_measurements.append(measurements)
    
    # Combine all measurements if we have any
    if all_measurements and save_plots:
        combined_measurements = pd.concat(all_measurements, ignore_index=True)
        combined_measurements.to_csv(f'{output_base_dir}/combined_measurements.csv', index=True)
        
        # Create summary statistics
        summary_stats = combined_measurements.groupby('image_name').agg({
            'area': ['count', 'sum', 'mean', 'std'],
            'axis_major_length': ['mean', 'std'],
            'axis_minor_length': ['mean', 'std'],
            'eccentricity': ['mean', 'std']
        }).round(3)
        
        summary_stats.to_csv(f'{output_base_dir}/summary_statistics.csv')
        
        print(f"\n--- BATCH PROCESSING COMPLETE ---")
        print(f"Results saved to: {output_base_dir}")
        print(f"Combined measurements: {len(combined_measurements)} total objects")
        print(f"Summary statistics saved for {len(image_files)} images")

# Example usage:

# For a single image:
# data = process_single_image('./path/to/image.jpg')

# For a folder of images:
process_folder('C:/Users/cirok/Desktop/Fix Segmentation/Input')

if __name__ == "__main__":
    # Example: Process a single image (replace with your image path)
    # data = process_single_image('./path/to/image.jpg')
    
    # Example: Process a folder of images (replace with your folder path)
    # process_folder('./my_image_folder')
    
    print("Script loaded with filled lumen functionality. Use process_single_image() or process_folder() functions.")