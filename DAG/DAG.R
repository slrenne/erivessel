library(pdftools)
library(magick)

# Define the PDF file path and the output PNG file path
pdf_file <- "DAG/DAG.pdf"
png_file <- "DAG/DAG.png"

# Convert the first page of the PDF to PNG (change density to adjust quality)
pdf_image <- image_read_pdf(pdf_file, density = 300)  # 300 DPI for high resolution
image_write(pdf_image, png_file, format = "png")
