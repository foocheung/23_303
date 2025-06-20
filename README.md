# In  Progress......



# Leveraging Optimized Oligonucleotide-Tagged Antigen Assemblies and Single-Cell Sequencing for Multiplexed Proteogenomic Profiling of Human B Cell Antigen Reactivities

This repository contains the core code and pipelines used in the forthcoming manuscript:

**"Leveraging optimized oligonucleotide-tagged antigen assemblies and single-cell sequencing for multiplexed proteogenomic profiling of human B cell antigen reactivities."**

Included methods span preprocessing, annotation, BCR clonotyping, CDR3 clustering, antigen specificity analysis, and cross-dataset clone comparison.

---

## Manuscript Reproducibility Repository

This repository provides all scripts and data needed to reproduce the figures and key analyses in the manuscript. Tools included support:

* BCR clonotyping and annotation
* CDR3 clustering and phylogenetic tree generation
* Clonal overlap analysis
* Antibody matching using the CoV-AbDab database

Note: Some figures may have undergone minor formatting adjustments for publication. These changes are purely aesthetic and do not affect the data.

---
## Folder Structure

The links below lead to folders that contain code to reproduce the corresponding figures and tables from the manuscript.

Repository structure with links to key folders:

- [ms_code/04.TableS1](https://github.com/foocheung/23_303/tree/main/ms_code/04.TableS1) – Scripts used to generate **Table S1**
- [ms_code/Figure3CD](https://github.com/foocheung/23_303/tree/main/ms_code/Figure3CD) – Code for **Figure 3C–D**
- [ms_code/Figure3E](https://github.com/foocheung/23_303/tree/main/ms_code/Figure3E) – Code for **Figure 3E**
- [ms_code/Figure5CDFG](https://github.com/foocheung/23_303/tree/main/ms_code/Figure5CDFG) – Code for **Figure 5C–D–F–G**
- [ms_code/Figure5E](https://github.com/foocheung/23_303/tree/main/ms_code/Figure5E) – Code for **Figure 5E**
- [ms_code/FigureS2](https://github.com/foocheung/23_303/tree/main/ms_code/FigureS2) – Code for **Supplementary Figure S2**
- [ms_code/FigureS4](https://github.com/foocheung/23_303/tree/main/ms_code/FigureS4) – Code for **Supplementary Figure S4**
- [ms_code/FigureS5](https://github.com/foocheung/23_303/tree/main/ms_code/FigureS5) – Code for **Supplementary Figure S5**
- [ms_code/FigureS6](https://github.com/foocheung/23_303/tree/main/ms_code/FigureS6) – Code for **Supplementary Figure S6**

Other folders:

- [data/](https://github.com/foocheung/23_303/tree/main/data) – Input files used for figure generation
- [README.md](https://github.com/foocheung/23_303/blob/main/README.md) – Overview of methods, tools, and usage instructions



---

## Repository Contents

| Path        | Description                                                                |
| ----------- | -------------------------------------------------------------------------- |
| `ms_code/`  | R scripts for BCR clonotyping, CDR3 clustering, and CoV-AbDab search       |
| `data/`     | Input files including clone tables, overlap matrices, and tree assignments |
| `README.md` | Overview of methods, tools, and repository usage                           |

---




## Quick Start

This guide will help you set up the repository and reproduce specific figures or tables from the manuscript. Each figure has its own script located in the `ms_code/` directory. To get started, clone the repository, navigate to the desired figure folder, install any required R libraries (as listed at the top of each script), and run the script to generate the corresponding outputs.

### 1. Clone the Repository

```bash
git clone https://github.com/foocheung/23_303
cd 23_303
````

### 2. Install Required R Packages

Open the R script for the figure you want to reproduce (e.g., `Figure3E`) and check the first few lines for `library()` calls. Install any missing packages using `install.packages()` for CRAN packages or `BiocManager::install()` for Bioconductor packages.

Example:

```r
install.packages(setdiff(c("readr", "dplyr", "ggplot2", "gridExtra", "ggrepel"), rownames(installed.packages())))
invisible(lapply(c("readr", "dplyr", "ggplot2", "gridExtra", "ggrepel"), library, character.only = TRUE))
```

> Tip: This approach ensures you only install the packages needed for each figure.

### 3. Run a Figure Script

Navigate into the folder for the figure you want to reproduce and run the script:

```bash
cd ms_code/Figure3E
Rscript covabdab-clone-scatter.R
```

Output files (plots, tables, summaries) will be saved in the same directory.

### 4. Input Data

All required input files are provided in the [`data/`](https://github.com/foocheung/23_303/tree/main/data) directory. No additional preprocessing is needed.



---




## Key Tools and Packages

* **Cell Ranger** (10x Genomics) – Alignment and quantification
  [https://www.10xgenomics.com/products/single-cell-gene-expression](https://www.10xgenomics.com/products/single-cell-gene-expression)

* **SingleR** – Reference-based cell type annotation
  [https://bioconductor.org/packages/SingleR/](https://bioconductor.org/packages/SingleR/)

* **CITE-seq-Count** – Antibody barcode quantification
  [https://github.com/Hoohm/CITE-seq-Count](https://github.com/Hoohm/CITE-seq-Count)

* **Immcantation v4.4.0** – BCR annotation and lineage analysis
  [https://immcantation.readthedocs.io/](https://immcantation.readthedocs.io/)

* **CoV-AbDab** – SARS-CoV-2 antibody database
  [https://opig.stats.ox.ac.uk/webapps/covabdab/](https://opig.stats.ox.ac.uk/webapps/covabdab/)

---

## CHI-Developed Tools and Pipelines

* **[upsetplotgen](https://github.com/foocheung/upsetplotgen)**  
  Golem-based Shiny app for BCR clonal overlap visualization and summary.  
  *Developed by Foo Cheung (2025)*

* **[fastaTreeClustering](https://github.com/foocheung/fastaTreeClustering)**  
  Interactive app for CDR3 clustering and tree visualization using Levenshtein distance.  
  *Developed by Foo Cheung (2025)*

* **[Demultiplexing_SNPSv2](https://github.com/foocheung/Demultiplexing_SNPSv2)**  
  Snakemake workflow for RNA-based SNP demultiplexing of pooled scRNA-seq data.  
  *Developed by Foo Cheung (2022)*

---



