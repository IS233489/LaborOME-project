# LaborOME-project


## Project Summary
![Labor-ome Study_Methods figure-022221-PILOT](https://raw.githubusercontent.com/IS233489/LaborOME-project/master/logos/Labor-omeStudy_Methodsfigure-022221-PILOT.jpg) Figure 1. Pilot project overview

We report on a cross-sectional investigation using a work-phase framework and multi-site composite skin sampling (Figure 1) to interrogate the skin microbiome, resistome, and mobilome of swine workers during a typical work shift that includes showering as a biosecurity intervention. We observed a significant increase in bacterial DNA load on worker skin during the workday shift, with concurrent changes in the composition and abundance of specific microbial taxa, ARGs, and mobile genetic elements (MGEs) that could harbor ARGs. We further observed that compulsory showering at the end of the workday reverted the skin microbiome and resistome to a baseline state, however the worker mobilome remained enriched with numerous key mobile genetic elements of public health significance. These results suggest that occupational work in swine facilities can significantly impact workers’ skin microbiomes, but that these impacts may be transient if biosecurity interventions are implemented. The relevance of these findings for short- and long-term swine worker health requires further study. However, the observation that showering may dampen daily microbiome impacts has important public health implications, as it demonstrates that biosecurity protocols could be leveraged to minimize microbial transmission from farms to general communities.

## Description of Analytical Scripts

# Microbiome Analysis of Swine Worker Study

This R script provides a comprehensive pipeline for integrating metadata, quality controlling 16S rRNA amplicon data, performing statistical modeling, and conducting downstream microbiome analyses including diversity estimation, contamination filtering, and relative abundance profiling.

---

## 📦 Required R Packages

This script depends on numerous CRAN and Bioconductor packages, including:

- `phyloseq`, `microbiome`, `microbiomeutilities`, `DESeq2`, `vegan`, `ggplot2`, `lme4`, `lmerTest`, `lsmeans`, `ALDEx2`, `decontam`, `SpiecEasi`, `MicrobiotaProcess`, `microViz`, `pairwiseAdonis`, `ggthemes`, `ggpubr`, `ggrepel`, `patchwork`, `cowplot`, `randomcoloR`, `htmltools`, and others.

Use `BiocManager::install()` and `install.packages()` as needed.

---

## 📁 Input Files

- `MCOHS_WORKER_METADATA.csv`: Initial participant metadata
- `Noyes_Project_032_MiSeq_Summary.csv`: Sequencing summary (e.g. read quality, counts)
- `track_asvtab_edited68samples.csv`: DADA2 tracking output
- `FINAL_MCOHS_WORKER_METADATA.csv`: Final merged metadata
- ASV table, taxonomy table, and metadata table (for phyloseq construction)
- `ilya_phyloseq_68samples_final_7.20.rds`: Initial phyloseq object

---

## 🧪 Key Analyses

### Metadata Integration
- Joins multiple metadata sources into a unified data frame
- Cleans numeric fields and prepares for modeling

### Exploratory Statistics
- Linear mixed models (LMMs) analyze effects of:
  - Collection Phase
  - Occupational Task
  - Task Exposure
- Outcomes analyzed: raw read counts, 16S qPCR values, read quality
- Visualized via `ggplot2` boxplots, jitter plots, and LSMEANS comparisons

### Phyloseq Construction
- Loads ASV, taxonomy, and sample data
- Constructs and saves a `phyloseq` object
- Filters out eukaryotic and USDA-blind samples
- Decontamination via `decontam` (qPCR-based frequency method)

### Microbiome Diversity & Relative Abundance
- Rarefaction curves via `MicrobiotaProcess`
- Taxonomic agglomeration and relative abundance plots
- Visualization of top phyla and families by `CollectionPhase`

### Network Analysis
- ASV filtering and subset by `CollectionPhase`
- Network estimation using `SpiecEasi`

---

## 💾 Output Files

- `*_plot.pptx`: Visualizations exported for publication
- `.rds` files: Saved `phyloseq` objects (filtered, grouped)
- `final_taxonomy_table.csv`: Exported taxonomy data

---

## 📌 Notes

- Use `set.seed()` when replicating rarefaction or network inference steps.
- Visualization uses color schemes from `ggthemes` Tableau palettes.

# Resistome and Mobilome Analysis Pipeline

This R script conducts quality control, statistical modeling, and taxonomic and functional profiling of shotgun metagenomic data focused on the resistome and mobilome from swine-exposed workers and environmental samples.

---

## 📦 Required Packages

- CRAN: `ggplot2`, `dplyr`, `ggpubr`, `ggrepel`, `patchwork`, `cowplot`, `lme4`, `lmerTest`, `lsmeans`, `scales`, `randomcoloR`, `htmltools`
- Bioconductor: `phyloseq`, `DESeq2`, `ALDEx2`, `microbiome`, `apeglm`, `ashr`, `metagenomeSeq`
- GitHub: `microbiomeutilities`, `MicrobiotaProcess`, `phylosmith`, `ranacapa`, `metagMisc`
- Others: `pairwiseAdonis`, `microViz`, `vegan`

Use `install.packages()` or `BiocManager::install()` as appropriate.

---

## 📁 Input Files

- `FINAL_SHOTGUN_MCOHS_WORKER_METADATA.csv`: Sample and metadata file
- Resistome matrix (AMR gene counts): CSV format with genes as rows and samples as columns
- Gene annotation taxonomy file
- Sample metadata for `phyloseq` construction

---

## 🔍 Key Steps

### 1. **Metadata Integration and Preprocessing**
- Loads and factors metadata (collection phases, tasks, exposure)
- Computes derived variables (e.g., host-removed reads)

### 2. **Read Depth and Quality Assessment**
- Visualizes sequencing depth pre/post trimming and host-removal
- Linear mixed models (LMM) assess differences by sample type (`CollectionPhase`, `OccupationalTask`)

### 3. **Phyloseq Object Construction**
- Builds `phyloseq` object with OTU, taxonomy, and sample data
- Subsets for medically important AMR genes

### 4. **Relative Abundance and Composition**
- Visualizes abundance by `Type`, `Class`, `Mechanism`, and `Group`
- Filters low-abundance groups
- Outputs stacked bar plots by `CollectionPhase` and `Sample`

### 5. **Alpha Diversity Analysis**
- Computes diversity metrics (`Observed`, `Shannon`, `Simpson`)
- LMMs evaluate diversity by metadata categories
- Outputs boxplots and statistics

### 6. **Prevalence and Abundance Analysis**
- Computes prevalence vs. abundance per ARG group
- Stratified by `CollectionPhase`
- Bubble plots and prevalence scatterplots

### 7. **Rarefaction Analysis**
- Uses `MicrobiotaProcess` for rarefaction curve plotting
- Calculates rarefaction-adjusted richness

---

## 📊 Outputs

- `.pptx` plots: Barplots, boxplots, rarefaction curves
- `.rds` files: Saved phyloseq objects
- `.csv`: Taxonomy tables, diversity metrics, prevalence matrices

---

## 💡 Notes

- Be sure to preprocess your AMR matrix to remove `RequiresSNPConfirmation` rows if applicable.
- Colors and themes use Tableau palettes for accessibility.
- Some parts of the script are modular; comment/uncomment sections as needed.


# Mobilome Analysis of Shotgun Metagenomes

This R script performs mobilome profiling using metagenomic data from human and environmental samples. It focuses on mobile genetic elements (MGEs) including plasmids, transposons, insertion sequences (IS), prophages, integrative conjugative elements (ICE), and viruses.

---

## 📦 Required R Packages

This script uses packages from CRAN, Bioconductor, and GitHub:

- **Core**: `phyloseq`, `ggplot2`, `dplyr`, `ggthemes`, `ggpubr`, `patchwork`, `vegan`, `lme4`, `lmerTest`, `lsmeans`
- **Microbiome tools**: `metagenomeSeq`, `ALDEx2`, `zCompositions`, `pairwiseAdonis`, `metagMisc`
- **Visualization**: `ggrepel`, `randomcoloR`, `scales`, `vioplot`, `cowplot`

Install missing dependencies using `install.packages()` or `BiocManager::install()`.

---

## 📁 Input Files

- `AMR_analytic_matrix_updatedID.csv`: Normalized MGE abundance matrix
- `WORKER_MOBILOME_ANNOTATION_FINAL.csv`: Gene annotation file with MGE types
- `FINAL_SHOTGUN_MCOHS_WORKER_METADATA.csv`: Sample metadata file

---

## 🧪 Key Analyses

### 1. **Phyloseq Object Construction**
- Converts matrices to `otu_table`, `tax_table`, and `sample_data`
- Builds a unified `phyloseq` object for downstream analysis
- Subsets to remove environmental samples where applicable

### 2. **Mobilome Composition & Relative Abundance**
- Taxonomic agglomeration by `Type` (e.g., plasmids, ICE, TE)
- Visualization of relative abundances across collection phases
- Barplots at type and gene family levels

### 3. **Alpha Diversity**
- Richness and Shannon diversity metrics
- Boxplots with jittered sample points
- Linear mixed-effects models test differences by `CollectionPhase`

### 4. **Beta Diversity & Ordination**
- Imputation with `zCompositions` for compositional data
- Euclidean distances via `rclr`-transformed counts
- PCoA ordination and centroid plots with ellipses
- ANOSIM and pairwise PERMANOVA tests (`pairwise.adonis`)

### 5. **Mobilome Subtype Analyses**
Includes independent analyses and plots for:
- **Plasmid structural genes** (e.g., conjugation, efflux, relaxase)
- **IS families** (e.g., IS1, IS3, IS200/IS605)
- **Transposons (TEs)**, **ICE**, **Prophages**, and **Viruses**

Each subtype includes:
- Relative abundance barplots (by sample and collection phase)
- Diversity analyses and statistical modeling
- Ordination plots and statistical group comparisons

---

## 📊 Outputs

- `.pptx` or `.png` plots: Barplots, ordination plots, diversity boxplots
- `.rds`: Saved `phyloseq` objects
- `.csv`: Richness statistics, ordination distances, group comparisons

---

## 💡 Notes

- Script contains separate analyses for imputed and non-imputed MGE data
- Uses color palettes like Tableau or manually defined ones for clarity
- Collection phases analyzed include:
  - `WORKDAY_START`, `WORKDAY_END`, `POST_SHOWER`, `SWINE`, `ENVIRONMENT`

# SpiecEasi Network Analysis

This R notebook performs microbial ecological network analysis using the SPIEC-EASI algorithm on multiple environmental and biological sample types.

## 🧬 Purpose

To infer and compare microbial association networks from different sample environments (e.g., worker start/end of shift, post-shower, swine, and environment) using SPIEC-EASI with Meinshausen-Bühlmann neighborhood selection.

## ✨ Features

- Loads `phyloseq` objects from `.rds` files for several sample types
- Runs `spiec.easi()` with 999 repetitions to infer robust networks
- Saves inferred networks for downstream comparison and visualization

## 📦 Requirements

Install the following R packages:

```r
install.packages("ggsci")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("phyloseq", "microbiome"))
install.packages(c("ggplot2", "dplyr", "ggrepel", "scales", "SpiecEasi"))
```

## 📂 Data Input

This script expects `.rds` files for each sample type which can be extracted from phyloseq objects as detailed in the Microbiome script:

- `SE_ps.noncontam.worker.WORKSTART.f.rds`
- `SE_ps.noncontam.worker.WORKEND.f.rds`
- `SE_ps.noncontam.worker.POSTSHOWER.f.rds`
- `SE_ps.noncontam.worker.SWINE.f.rds`
- `SE_ps.noncontam.worker.ENVIRONMENT.f.rds`

## 💾 Output

Each SPIEC-EASI model is saved as an `.rds` file, such as:

- `SE_WORKSTART.mb999.rds`
- `SE_WORKEND.mb.rds`

## 🚀 Example Usage

```r
set.seed(1244)
SE_WORKSTART.mb <- spiec.easi(SE_ps.noncontam.worker.WORKSTART,
                              method='mb',
                              lambda.min.ratio=1e-1,
                              nlambda=20,
                              icov.select.params=list(rep.num=999))
saveRDS(SE_WORKSTART.mb, file = "SE_WORKSTART.mb999.rds")





# MAG Analysis and Visualization

This repository contains scripts for analyzing and visualizing Metagenome-Assembled Genomes (MAGs) based on quality metrics and metadata annotations.

## Overview

The `MAG_ANALYSIS_script.ipynb` performs the following:

- Activates a Conda environment and installs required R packages.
- Loads MAG metadata from a CSV file (`itol_tree_sorted_metadata_R_data.csv`).
- Converts selected metadata columns to numeric types for analysis.
- Creates annotated scatterplots of genome completeness versus contamination:
  - Points colored by `Collection phase`.
  - Points sized by `GC content`.
  - Quality thresholds highlighted with shaded regions.
  - Marginal density plots added for better distribution visualization.
- Saves the final figure as a high-resolution `.svg` file.

## Requirements

- Conda
- R environment
- Installed R packages:
  - `tidyverse`
  - `ggplot2`
  - `ggExtra`
  - `cowplot`
  - `svglite`
  - `ggpubr`
  - `dplyr`

You can install the R packages via:

```bash
conda activate r-MAG_ANALYSIS
conda install -c conda-forge r-tidyverse
```
(Other packages are typically included in `tidyverse`, or can be installed separately using `install.packages()`).

## Usage

1. Clone this repository:
    ```bash
    git clone https://github.com/yourusername/your-repo-name.git
    cd your-repo-name
    ```

2. Activate the Conda environment:
    ```bash
    conda activate r-MAG_ANALYSIS
    ```

3. Open the notebook:
    ```bash
    jupyter notebook MAG_ANALYSIS_script.ipynb
    ```

4. Ensure the input file `itol_tree_sorted_metadata_R_data.csv` is present in the same directory.

5. Run the notebook cells sequentially to perform the analysis.

6. The final plot will be saved as `COMP_CONT.svg` in the working directory.

## Output

- **COMP_CONT.svg**: A high-quality, annotated scatterplot with marginal distributions.

## Notes

- The script uses manual color mapping for collection phases.
- Quality thresholds (like >90% completeness and <5% contamination) are visually highlighted.
- Additional exploratory plots (draft genome size, species novelty) are partially included but may require further editing.

