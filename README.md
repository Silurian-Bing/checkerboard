# Methodological Approach for Spatial Distribution Analysis
Fossil distribution
This archive contains four files, each corresponding to a different step in the data processing and analysis workflow for spatial fossil distribution.

File Overview
Image_for_data.png
This image file contains the original data source. All specimen locations were carefully traced by hand from the physical specimen onto the screen, creating an accurate distribution map of the fossils.

Data_Acquisition.html
This file is an image recognition tool written in JavaScript and HTML. Open this file in a web browser to digitize the distribution points from Image_for_data.png. The program accurately converts the traced locations into pixel coordinates and, most importantly, automatically measures the pairwise distances between all specimens in the illustration. This automated measurement replaces what would otherwise be a prohibitively time-consuming manual task. The resulting data are saved as a CSV file.

Dataset_S1.csv
This CSV file contains the processed spatial data exported from the Data_Acquisition.html program. It includes the pixel coordinates and pairwise distances for all specimens.

spatial_analysis.R
This R script is used for further analysis of the dataset, including the calculation of Voronoi (Thiessen) polygons and other spatial statistics. The results presented in the main paper and Supplementary Information (SI) were obtained using this script.

Workflow Summary
Manual digitization: Specimen positions are hand-traced onto Image_for_data.png.
Automated data extraction: Use Data_Acquisition.html to digitize positions and automatically calculate all pairwise distances, saving the results as Dataset_S1.csv.
Statistical analysis: Analyze the dataset using spatial_analysis.R to generate Voronoi polygons and final results.
For further details, refer to the main paper and Supplementary Information.
