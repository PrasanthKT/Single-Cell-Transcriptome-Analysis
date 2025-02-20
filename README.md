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

### Getting Started
1.	Install all the dependencies.
2.	Download the count matrix files from GEO
3.	Navigate to the Scripts folder and run the R script
4.	Use the output files from the analysis for further processing and analysis in R.

### Workflow
The detailed workflow for the analysis is documented in the report, which is available in the repository as ```Report.pdf```. Shell scripts for each step of the workflow are organized under the ```Scripts``` folder, while R scripts and their corresponding input files are located in the ```R_Analysis```folder within the repository.

### Results
A comprehensive interpretation of the results is provided in the ```Report.pdf``` report. Key findings are summarized in the rport.
