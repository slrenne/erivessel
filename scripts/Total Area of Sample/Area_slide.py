import cv2
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
from tqdm import tqdm 

def calculate_tissue_area(image_path, save_visualizations=False, output_folder=None):
    # Read the image
    img = cv2.imread(image_path)
    if img is None:
        print(f"Error: Could not read image {image_path}")
        return 0, None, None
    
    # Convert to RGB (OpenCV loads as BGR)
    img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    
    # Convert to LAB color space
    lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
    
    # Extract the L channel (luminance)
    l_channel = lab[:,:,0]
    
    # Apply threshold to separate tissue from background
    _, binary = cv2.threshold(l_channel, 220, 255, cv2.THRESH_BINARY_INV)
    
    # Clean up the binary image
    kernel = np.ones((5,5), np.uint8)
    binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel)
    binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, kernel)
    
    # Find contours
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Create a mask for the full tissue
    mask = np.zeros_like(l_channel)
    
    # Draw all contours above a certain size (to avoid small artifacts)
    min_contour_area = 1000  # Adjust based on your image size
    for contour in contours:
        if cv2.contourArea(contour) > min_contour_area:
            cv2.drawContours(mask, [contour], 0, 255, -1)
    
    # Calculate area in pixels
    area_pixels = cv2.countNonZero(mask)
    
    # Create visualization of the filled tissue
    filled_tissue = cv2.bitwise_and(img_rgb, img_rgb, mask=mask)
    
    if save_visualizations:
        if output_folder is None:
            output_folder = os.path.dirname(image_path) or '.'
        
        # Create visualization
        plt.figure(figsize=(15, 10))
        
        plt.subplot(2, 2, 1)
        plt.imshow(img_rgb)
        plt.title('Original Image')
        plt.axis('off')
        
        plt.subplot(2, 2, 2)
        plt.imshow(binary, cmap='gray')
        plt.title('Thresholded Image')
        plt.axis('off')
        
        plt.subplot(2, 2, 3)
        plt.imshow(mask, cmap='gray')
        plt.title('Tissue Mask')
        plt.axis('off')
        
        plt.subplot(2, 2, 4)
        plt.imshow(filled_tissue)
        plt.title(f'Filled Tissue (Area: {area_pixels} pixels)')
        plt.axis('off')
        
        plt.tight_layout()
        
        # Get filename without extension
        base_name = os.path.splitext(os.path.basename(image_path))[0]
        plt.savefig(f"{output_folder}/{base_name}_analysis.png", dpi=300)
        plt.close()
    
    return area_pixels, mask, filled_tissue

def process_folder(input_folder, output_folder, file_extensions=None):
    if file_extensions is None:
        file_extensions = ['.jpg', '.jpeg', '.png', '.tif', '.tiff']
    
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Get all image files in the input folder
    image_files = []
    for ext in file_extensions:
        image_files.extend(glob.glob(os.path.join(input_folder, f'*{ext}')))
        image_files.extend(glob.glob(os.path.join(input_folder, f'*{ext.upper()}')))
    
    if not image_files:
        print(f"No image files found in {input_folder} with extensions {file_extensions}")
        return None
    
    print(f"Found {len(image_files)} image files to process")
    
    # Process each image and collect results
    results = []
    
    for img_path in tqdm(image_files, desc="Processing slides"):
        try:
            # Get filename without path
            filename = os.path.basename(img_path)
            
            # Calculate area
            area, mask, filled_tissue = calculate_tissue_area(img_path)
            
            # Save the original image
            img = cv2.imread(img_path)
            base_name = os.path.splitext(filename)[0]
            cv2.imwrite(f"{output_folder}/{base_name}_original.png", img)
            
            # Save the mask
            if mask is not None:
                cv2.imwrite(f"{output_folder}/{base_name}_mask.png", mask)
            
            # Store the results
            results.append({
                'Image': filename,
                'Tissue Area (pixels)': area
            })
            
        except Exception as e:
            print(f"Error processing {img_path}: {e}")
    
    # Create a DataFrame and save to CSV
    if results:
        df = pd.DataFrame(results)
        csv_path = os.path.join(output_folder, "tissue_areas.csv")
        df.to_csv(csv_path, index=False)
        print(f"Results saved to {csv_path}")
        return df
    else:
        print("No results were generated")
        return None

# Example usage
if __name__ == "__main__":

    # Process the arguments
    input_folder = 'Input'
    output_folder = 'Output'
    
    # Run the analysis
    print(f"Processing slides from {input_folder}")
    print(f"Results will be saved to {output_folder}")
    
    df = process_folder(input_folder, output_folder)
    
    if df is not None:
        print("Processing complete!")
        print(f"CSV file with area measurements saved to {output_folder}/tissue_areas.csv")
    else:
        print("No results were generated.")