# Morphological-traits-and-the-frugivory-degree-in-avian-fruit-consumers
# Morphological traits and the frugivory degree in avian fruit consumers

## Overview

This repository contains the data and R code supporting a study investigating the relationships between morphological traits and the degree of frugivory in birds. The analyses examine variation in body size, beak morphology, and digestive organ traits, as well as their associations with frugivory and fruit consumption patterns.

The repository also includes phylogenetic analyses to account for evolutionary relationships among bird species and to investigate phylogenetic patterns in plant–bird interactions.

## Repository contents

| File | Description |
|---|---|
| `Crestani_Biotropica_GitHub_updated.R` | R script containing the statistical and phylogenetic analyses |
| `morpho_dados_25.csv` | Individual-level morphological measurements of birds |
| `morpho_rel_dados_25.csv` | Morphological and frugivory data used in the statistical analyses |
| `interactions_25.csv` | Plant–bird interaction matrix containing fruit consumption rates |
| `bird_phy_v05.nex` | Bird phylogenetic tree |
| `crestani_08_2025.nex` | Set of bird phylogenetic trees |
| `plant_tree.nex` | Plant phylogenetic tree |

## Statistical analyses

The R script includes the following analyses:

- Correlations among morphological traits.
- Allometric relationships between body mass and morphological traits.
- Beta regression examining the relationship between morphological traits and the degree of frugivory.
- Principal component analyses (PCA) of morphological traits and fruit consumption patterns.
- Procrustes analysis comparing morphological and interaction ordinations.
- Phylogenetic signal analyses of bird morphological traits and frugivory.
- Phylogenetic generalized least squares (PGLS) models.
- Phylogenetic signal analysis of plant–bird interactions.

## Data and code availability

The data and R code are publicly available in this repository:

https://github.com/anacristinacrestani/Morphological-traits-and-the-frugivory-degree-in-avian-fruit-consumers

## Reproducibility

To reproduce the analyses, download or clone this repository, install the R packages required by the script, and run `Crestani_Biotropica_GitHub_updated.R` in R or RStudio.

The input files should be kept in the same working directory as the script, unless the file paths are adjusted accordingly.

## Citation

If you use these data or analyses, please cite the associated research article.

Full publication details will be added when available.
