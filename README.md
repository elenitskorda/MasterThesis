# Master Thesis Title

# Finding subtypes of ovarian cancer in omics data for precision medicine

### Name and e-mail addresses

### **Student:** Eleni Theofania Skorda, <el5660sk-s@student.lu.se>

### **Supervisor:** Fredrik Levander, <fredrik.levander@immun.lth.se>

## Specific Aim

 The specific aim of this thesis is the identification of molecular subtypes of high-grade ovarian cancer using
 proteomics, phosphoproteomics, and other omics data and to explore potential personalized treatments for
 each subtype.

### Abbrevations

## Methods

* **Data collection:** Identify and obtain datasets for CPTAC data portal that include proteomics, phosphoproteomics, and other omics data for high-grade ovarian cancer. Additionally local cell line data are available
 for comparison.
* **Data preprocessing:** Clean and preprocess the data, normalize the data, and address batch effects.
* **Unsupervised clustering:** Use unsupervised clustering methods to identify potential subtypes of highgrade ovarian cancer based on the omics data.
* **Differential abundance analysis:** Perform differential abundance analysis to identify the proteins or
 phosphoproteins that are significantly differentially expressed between the identified subtypes of high-grade
 ovarian cancer.
* **Integrating cell line data:** Explore strategies to integrate the omics data with local cell line data, both
 at the pathway and individual levels.
* **Signature development:** Use machine learning algorithms to develop signatures that distinguish between
 the identified subtypes of high-grade ovarian cancer. Downstream analysis with python
* **Multi-omics data exploration:** Integrate multiple omics data types (e.g., proteomics, phosphoproteomics,
 transcriptomics, etc.) to identify potential targets for precision medicine.
* **Spatial omics data:** Explore the spatial omics data to identify spatially heterogeneous subpopulations of
 cells that may contribute to the identified subtypes of high-grade ovarian cancer

## Requirements

* Fragpipe - version 20.0
* R programming language - version 4.3.1
* OmicLoupe - version 1.1.7

### R Packages/Libraries
to be continued ...
### Input files

> The Data are coming from National Cancer Institute's s Office of Cancer Clinical Proteomics Research (OCCPR). 
More specifically, data are coming from Proteomic Data Commos (PDC) which represents the public repository of proteogenomic
tumor datasets.  Proteomic data and related data files are organized into datasets by tumor type, study and sub-proteome. 
In addition to the raw mass spectrometry-based data files, computational processing is performed to map spectra to peptide 
sequences and identify proteins. All spectra collected in CPTAC studies are analyzed for consistency and reproducibility as 
part of the CPTAC Common Data Analysis pipeline and an interactive QC report for each study is generated. All data is freely 
available to the public.

### Computational Power

Doris + mine