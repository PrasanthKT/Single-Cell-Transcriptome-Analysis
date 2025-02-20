# Single-Cell-Transcriptome-Analysis

### Objective
The study aims to reproduce the original study from the paper ("Single-cell RNA-seq reveals a resolving immune phenotype in the oral mucosa") to use single-cell RNA sequencing (scRNA-seq) to characterize immune cell heterogeneity in the oral mucosa, identifying distinct subsets and their gene expression profiles. It focuses on understanding immune responses in inflammation, such as periodontitis, and their resolution post-treatment

### Tools and Dependencies
1. Seurat
2. R and Bio conductor Packages

### Installation

R and Bioconductor packages, install them directly in R:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("Seurat"))
```

### Input files 

GEO ID:GSE244633
File Type: Count matrix files (HDF5 format)

### Getting Started
1.	Install all the dependencies.
2.	Download the count matrix files from GEO
3.	Navigate to the Scripts folder and run the R script

### Results
A comprehensive interpretation of the results is provided in the ```Report.pdf``` report. Key findings are summarized in the rport.
