# Navigating floristic networks: linking oceanic currents and littoral flora across oceanic archipelagos

This repository contains all files and analyses needed to reproduce the results of our study. 

## Folder Structure

The repository follows the following folder structure:

```
- a_script/: Contains the source code files for our analyses.
  /Functions/: Contains the stand-alone functions we call in our scripts.
- b_data/: Contains the data files for the project.
- c_outputs/: Contains all the outputs of our analyses, including figures.
- image_editor/: contains .svg files for editing certain figures.
- Qgis/: contains Qgis files to retrieve specific data.
```

The analyses require to be run in a certain sequence:

1) **compute_floristic_distances.Rmd**: reads floristic data from databases and computes floristic distance matrices based on Bray Curtis DI.

2) **compute_current_connectivity.Rmd**: downloads HYCOM data and uses it to compute pairwise cost distances by currents between islands in each archipelago.

> The script downloads a large amount of NetCDF files. An option is to download these directly from the [Zenodo] folder and read them locally: 
https://zenodo.org/records/12659307?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjgxMzgwODk4LTQ5NDktNDgwMC1iZTYwLWIyNDE2NjY5NTQ4NSIsImRhdGEiOnt9LCJyYW5kb20iOiIyNGU5MDg3NTAyM2E4MmZhNTllZmRjMjdlNDAyZWYyMyJ9.dgtOQZgvOw5vSmq-hyK6uJteSIeAL80YHtCYs69z0HSKAbL5_45m_sZpyRuCaQPxuDRmBzCujqppeXfI6E57lw


3) **correlations.Rmd**: reads floristic and current distance matrices produced by scripts in (1) and (2) and computes Procrustes correlations between them.

4) **network_analyses.Rmd**: plots raw data on % littoral plants by islands, computes centrality measures of each island based on their position in the current connectivity network, runs GLM models to test whether % littoral plants is predicted by centrality measures, and plots networks.

5) **plot_single_day_connectivity.Rmd**: plots an example of least cost path between two points in the Galapagos archipelago based on current connectivity.