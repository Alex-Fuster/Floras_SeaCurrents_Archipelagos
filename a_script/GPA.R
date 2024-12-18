
library(vegan)
library(FactoMineR)
library(ade4)
library(vegan)
library(tidyverse)
library(cluster)
library(dplyr)
library(tidyverse)
library(cluster)
library(ggplot2)
library(reshape2)
library(viridis)
library(ggpubr)
source(here::here("a_script/Functions.R"))


######################################### Obtain data

# Load floristic distances

dist_mat_all_litt_gal <- as.matrix(read_rds(here::here("c_output/floristic_dist_matrices/dist_mat_all_litt_gal.rds")))

#Order rows and columns and convert matrices to distance objects

dist_mat_all_litt_gal <- dist_mat_all_litt_gal[order(rownames(dist_mat_all_litt_gal)), order(colnames(dist_mat_all_litt_gal))]

# Load current distances
mat_medianMin_sym_Galapagos <- read_rds(here::here("c_output/current_conn_matrices/mat_medianMin_sym_Galapagos.rds"))


# Load currents geographic matrices and convert to distance objects

mat_geo_gal <- read.csv(here::here("b_data/Data_Galapagos/geographic_distance_galap.csv"))
rownames(mat_geo_gal) <- mat_geo_gal[,1]
mat_geo_gal<-mat_geo_gal[,-1]

## order names and convert to matrices

mat_geo_gal <- as.matrix(mat_geo_gal[order(rownames(mat_geo_gal)), order(colnames(mat_geo_gal))])



#####################################################


# Perform NMDS on the floristic distance matrix
nmds_floristic <- metaMDS(dist_mat_all_litt_gal, k = 2, distance = "euclidean")  
floristic_axes <- as.data.frame(scores(nmds_floristic))  # Extract NMDS axes

# Perform NMDS on the geographic distance matrix
nmds_geo <- metaMDS(mat_geo_gal, k = 2, distance = "euclidean")  
geo_axes <- as.data.frame(scores(nmds_geo))  # Extract NMDS axes

# Perform NMDS on the currents distance matrix
nmds_currents <- metaMDS(mat_medianMin_sym_Galapagos, k = 2, distance = "euclidean")  
currents_axes <- as.data.frame(scores(nmds_currents))  # Extract NMDS axes

# Optional: Check stress values to evaluate the fit
nmds_floristic$stress
nmds_geo$stress
nmds_currents$stress


# Combine NMDS axes into a single data frame
combined_axes <- data.frame(
  floristic_axis1 = floristic_axes$NMDS1,
  floristic_axis2 = floristic_axes$NMDS2,
  geo_axis1 = geo_axes$NMDS1,
  geo_axis2 = geo_axes$NMDS2,
  currents_axis1 = currents_axes$NMDS1,
  currents_axis2 = currents_axes$NMDS2
)

rownames(combined_axes) <- rownames(dist_mat_all_litt_gal)

# Perform Generalized Procrustes Analysis (GPA)
gpa_result <- GPA(combined_axes[-1,], group = c(2, 2, 2))  # Grouping the axes (2 axes per group)

# Print the GPA result
print(gpa_result)

# Check Procrustes Similarity Indexes
gpa_result$simi

# Check RV Coefficients
gpa_result$RV

# Check Procrustes ANOVA results
gpa_result$PANOVA


plot(gpa_result)


##################################

# function to perform NMDS with different dimensions
nmds_stress_plot <- function(dissimilarity_matrix) {
  max_dimensions <- 6  # Set the maximum number of dimensions to test
  stress_values <- numeric(max_dimensions)
  
  for (k in 1:max_dimensions) {
    nmds_result <- metaMDS(dissimilarity_matrix, k = k, trymax = 100)
    stress_values[k] <- nmds_result$stress
  }
  
  # Plot stress vs number of dimensions
  plot(1:max_dimensions, stress_values, type = "b", 
       xlab = "Number of Dimensions", ylab = "Stress Value", 
       main = "NMDS Stress vs Number of Dimensions")
  abline(h = 0.1, col = "red", lty = 2)  # Add a reference line at stress = 0.1
}

nmds_stress_plot(dist_mat_all_litt_gal)
nmds_stress_plot(mat_geo_gal)
nmds_stress_plot(mat_medianMin_sym_Galapagos)
