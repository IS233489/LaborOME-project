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

- Many installation lines are included in the script for reproducibility; you may want to comment them out once packages are installed.
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




