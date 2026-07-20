# IB3D
Three-Dimensional Ecological Connectivity Model for Landscape Genetics


This repository contains the R-based computational framework for **IB3D**, a 3D resistance model designed for landscape genetics analysis in marine environments. 

The pipeline computes effective resistance ($R_{eff}$) across a 3D grid using environmental variables (bathymetry and noise-attenuated shipping data) to model spatial connectivity and genetic isolation in marine megafauna populations.

## Key Features
* **3D Grid Architecture:** Incorporates vertical resolution (5m layers) alongside spatial coordinates.
* **Circuit Theory Approach:** Utilizes Kirchhoff’s Current Law to calculate effective resistance matrix.
* **Memory-Optimized Solver:** Replaces dense Cholesky decomposition with a scalable Conjugate Gradient (CG) algorithm for large pairwise datasets.
* **GIS Integration:** Automated pipeline to export localized isolation indices into standardized GeoTIFF layers via `terra`.

## Requirements
* R (>= 4.0)
* Required packages: `terra`, `Matrix`

## Usage
1. Place your environmental rasters (`bathymetry.tif`, etc.) and sample coordinates in the working directory.
2. Run the core script to generate the $63 \times 63$ pairwise resistance matrix.
3. Export the spatial connectivity layers for visualization in QGIS.

## Citation
If you use this model or code in your research, please cite it as:
> *Martin Žirovnický (2026). IB3D: Three-Dimensional Ecological Connectivity Model for Landscape Genetics. GitHub repository.*
