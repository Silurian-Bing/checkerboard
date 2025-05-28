# Methodological Approach for Spatial Distribution Analysis
This archive contains five files, each corresponding to a different step in the data processing and analysis workflow for spatial fossil distribution.

File Overview
Image_for_data.png
This image file contains the original data source. All specimen locations were carefully traced by hand from the physical specimen onto the screen, creating an accurate distribution map of the fossils.

Data_Acquisition.html
This file is an image recognition tool written in JavaScript and HTML. Open this file in a web browser to digitize the distribution points from Image_for_data.png. The program accurately converts the traced locations into pixel coordinates and, most importantly, automatically measures the pairwise distances between all specimens in the illustration. This automated measurement replaces what would otherwise be a prohibitively time-consuming manual task. The resulting data are saved as a CSV file.

Tool.html
This file is an enhanced version of Data_Acquisition.html. In addition to all the features of Data_Acquisition.html, Tool.html provides integrated spatial statistical analysis. If you encounter issues running the R script (for example, due to R version incompatibility or package installation problems), you can use Tool.html as an alternative. By simply uploading a distribution map—such as Image_for_data.png—directly into the tool, you can obtain analysis results (e.g., Voronoi polygons, nearest-neighbor statistics) within your web browser, without the need for additional software installations.

Dataset_S1.csv
This CSV file contains the processed spatial data exported from the Data_Acquisition.html (or Tool.html) program. It includes the pixel coordinates and pairwise distances for all specimens.

spatial_analysis.R
This R script is used for further analysis of the dataset, including the calculation of Voronoi (Thiessen) polygons and other spatial statistics. The results presented in the main paper and Supplementary Information (SI) were obtained using this script.

Workflow Summary
Manual digitization:
Specimen positions are hand-traced onto Image_for_data.png.

Automated data extraction:
Use Data_Acquisition.html to digitize positions and automatically calculate all pairwise distances, saving the results as Dataset_S1.csv.

Statistical analysis:
Analyze the dataset using spatial_analysis.R to generate Voronoi polygons and final results.

Alternative integrated analysis:
If R analysis cannot be performed (e.g., due to version or package issues), open Tool.html in your web browser. Upload your distribution map (such as Image_for_data.png) to directly obtain spatial analysis results. This provides a convenient, platform-independent alternative for data analysis.
