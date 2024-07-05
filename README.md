---
# Finding subtypes of ovarian cancer in omics data for precision medicine

---
*README for Master Thesis in Bioinformatics (60 cp).*

**Student:** Eleni Theofania Skorda, <el5660sk-s@student.lu.se>
**Supervisor:** Fredrik Levander, <fredrik.levander@immun.lth.se>

## Specific Aim

The specific aim of this thesis is the identification of molecular subtypes of high-grade ovarian cancer using proteomics and phosphoproteomics to explore potential personalized treatments for each subtype. This README introduces all the scripts and data necessary to reproduce the analysis.

## Methods

* **Data collection:** Identify and obtain datasets for CPTAC data portal that include proteomics, phosphoproteomics, and other omics data for high-grade ovarian cancer. Additionally local cell line data are available for comparison.
* **Data preprocessing:** Clean and preprocess the data, normalize the data, and address batch effects.
* **Unsupervised clustering:** Use unsupervised clustering methods to identify potential subtypes of highgrade ovarian cancer based on the omics data.
* **Integrating cell line data:** Explore strategies to integrate the omics data with local cell line data, both at the pathway and individual levels.
* **Differential abundance analysis:** Perform differential abundance analysis to identify the proteins or phosphoproteins that are significantly differentially expressed between the identified subtypes of high-grade ovarian cancer.
* **Multi-omics data exploration:** Integrate multiple omics data types (e.g., proteomics, phosphoproteomics, transcriptomics, etc.) to identify potential targets for precision medicine.
* **Enrichment Analysis:** Perform enrichment analysis to identify common pathways between benign and no benign cell lines subsequent the integration of proteome that may contribute to the identified unique pathways of high-grade ovarian cancer.

## Requirements

* Fragpipe (version 20.0)
* R programming language (version 4.3.1)
* OmicLoupe (version 1.1.7)

## R Packages/Libraries
Some of the libraries that have been used in this thesis have been imported from Bioncoductor and some other from R as well as GitHub. You can find all of them below, as well as their versions:


### Installation from GitHub
The Seurat package used for the integration and OmicLoupe has to be installed over github as follows.
  ```{R}
  remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
  library(Seurat)
  setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
  install.packages("presto")
  
  devtools::install_github("ComputationalProteomics/OmicLoupe")
  ```

### Installed from Bioconductor

*  *VennDetail* (1.16.0)
*  *ComplexHeatmap* (2.16.0)
*  *limma* (3.57.4)
*  *PCAtools* (2.12.0)
*  *NormalyzerDE* (1.18.0)

 GitHub

  Installations were made in R with for example:
  ```
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("VennDetail")
  
  ```

### R Based
* *devtools* (2.4.5)
* *SeuratObject* (5.0.1)
* *ggtext* (0.1.2)
* *tidyverse* (2.0.0) 
* *writexl* (1.4.2)
* *enrichR* (3.2)
* *patchwork* (1.2.0)
* *metap* (1.9)
* *VennDiagram* (1.7.3)
* *ggtext* (0.1.2)
* *dplyr* (1.1.3)
* *matrixStats* (1.0.0)
* *RColorBrewer* (1.1-3)
* *pheatmap* (1.0.12)
* *factoextra* (1.0.7)
* *ggpubr* (0.6.0)


## Input files

The Data is coming from National Cancer Institute's s Office of Cancer Clinical Proteomics Research (OCCPR). More specifically, the data is coming from Proteomic Data Commos (PDC) which represents the public repository of proteogenomic tumor datasets.  Proteomic data and related data files are organized into datasets by tumor type, study and sub-proteome. 
In addition to the raw mass spectrometry-based data files, computational processing is performed to map spectra to peptide sequences and identify proteins. All spectra collected in CPTAC studies are analyzed for consistency and reproducibility as part of the CPTAC Common Data Analysis pipeline and an interactive QC report for each study is generated. All data is freely available to the public.The second source of data for this study is cell line data acquired from different banks. These are the American Type Culture Collection (ATCC), Sigma-Aldrich, Japanese Collection of Research Bioresources (JCRB) Cell Bank, Riken Cell Bank and Henrik Olsen's Lab. These cell lines are: CAOV3, NIHOVCAR3, ES2, SKOV3, A2780, IGROV1, OVSAHO, KURAMOCHI, JHOC5 (Control), HeLa (Control). Each cell line contains 6 samples.

### Local Cell Line data files

* *"OvCa Proteins 5sub Loess DE_stats.tsv"*: contains proteome of cell line data.
* *"Phospho.pg_matrix(protein level).tsv"*: contains phosphoproteome of cell line data.
* *"OvCa phosphoProteins-DE_Normalised.tsv.phosphosites.tsv"*: contains phosphoproteome as well as their phosphosite of cell line data. This set of data is already normalized whereas the others are not.

### CPTAC Data

* *"ratio_peptide_MD.tsv"*: contains phosphoproteome on peptide level of CPTAC data. (Normalized)
* *"abundance_gene_MD_18.tsv"*: contains proteome of CPTAC data. (Normalized)
* *"abundance_protein_MD.tsv"*: contains phosphoproteome on protein level of CPTAC data. (Normalized)
* *"abundance_single-site_MD.tsv"*: contains phosphoproteome on single site level of CPTAC data. (Normalized)

#### For Fragpipe:

Based on the experiment type and analytical fraction, the TMT/iTRAQ quantification tutorial has been used for multiple plexes with a pooled reference sample. The workflows used are TMT 10-bridge and TMT 10-phospho-bridge for proteome and phosphoproteome respectively.

#### For Integration:

* *[abundance/ratio]_gene\_[normalization].tsv* which contains isobaric quantification information summarized from the psm.tsv tables by TMT-Integrator to the gene level.
* *[abundance/ratio]\_protein\_[normalization].tsv* which contains isobaric quantification information summarized from the psm.tsv tables by TMT-Integrator to the protein level.
* *[abundance/ratio]\_single-site\_[normalization].tsv* which contains isobaric quantification information summarized from the psm.tsv tables by TMT-Integrator based on modification sites that have been observed and quantified together.
* Protein expression evaluated across 3 cell lines. See paper [*Integrative proteomic profiling of ovarian cancer cell lines reveals precursor cell associated proteins and functional status F. Coscia, K. M. Watters, M. Curtis, M. A. Eckert, C. Y. Chiang, S. Tyanova, A. Montag, R. R. Lastra, E. Lengyel & M. Mann*](https://www.nature.com/articles/ncomms12645).

## Analysis

1. NormalyzerDE: *Stats.R*
2. OmicLoupe: *Stats.R*
3. Pre-processing Analysis: *unsupervised_clustering.R*
4. Unsupervised clustering: *unsupervised_clustering.R*
5. Integration Analysis: *proteome_integration.R*, *phosphoproteome_integration_protein.R*, *phosphoproteome_integration_singlesite.R*, *validation_of_integration.R*
7. Pathway analysis: *proteome_integration.R*
8. Differential Gene Expression Analysis: *proteome_integration.R*

Note that filepaths have to be adjusted in the scripts if the analysis is to be replicated. Further all input files have to be in the same filepath as the script to be executed.
