# -*- coding: utf-8 -*-

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Date: Mon Jun 5 10:28:18 2024                                                      ::
##  File Name: validation.R                                                            ::
##  Author: Eleni Theofania Skorda                                                     ::
##                                                                                     ::
##  Description:                                                                       ::
##      Integrates proteomic local cell line data (Reference) with CPTAC data (Query)  ::
##      in order to validate lines within the query dataset.                           ::
##                                                                                     ::
##  List of Functions:                                                                 ::
##      none                                                                           ::
##                                                                                     ::
##  Procedure:                                                                         ::
##      1. Prepare reference and query datasets to a proper format                     ::
##      2. Performing integration                                                      ::
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library("writexl")
library(tidyverse)
library(patchwork)

##==========================================================================
##  SECTION 1: Prepare reference and query datasets to a proper format    ==
##==========================================================================


#                           QUERY DATASET
##---------------------------------------------------------------------------

# Load the data into R
data <- read.delim("reference_paper.txt", sep = "\t", check.names = FALSE)

# Extract gene names and transpose the data
Gene = data[, 1]
rownames(data) = make.names(Gene, unique = TRUE)   
data = data[, -1]
data = as.data.frame(data)

# Transform data back to Raw values
data=2^data

## Impute low value for missing data in query
# Identify numeric columns
numeric_cols <- sapply(data, is.numeric)

# Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(data[, i], na.rm = TRUE)
  
  # Define a low value 
  low_val <- min_val
  
  #  # Replace NA values with the low value
  data[is.na(data[, i]), i] <- low_val
}

##---------------------------------------------------------------------------
#                              REFERENCE DATASET
##---------------------------------------------------------------------------

# Read refence data and remove unnecessary columns
local.cell.line <- read.delim("OvCa Proteins 5sub Loess DE_stats.tsv",
                      sep = "\t", check.names = FALSE)

# Remove unnecessary columns from the query
local.cell.line = subset(local.cell.line, select = -c(1:3))
local.cell.line = subset(local.cell.line, select = -c(2:12))

# Extract gene names and transpose the query
Gene = local.cell.line[, 1]
rownames(local.cell.line) = make.names(Gene, unique = TRUE)
local.cell.line = local.cell.line[, -1]

# Transform data back to Raw values
local.cell.line=2^local.cell.line

## Impute low value for missing data in reference
# Identify numeric columns
numeric_cols <- sapply(local.cell.line, is.numeric)

# Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(local.cell.line[, i], na.rm = TRUE)
  
  # Define a low value (e.g., min_val - 1)
  low_val <- min_val
  
  #  # Replace NA values with the low value
  local.cell.line[is.na(local.cell.line[, i]), i] <- low_val
}

##---------------------------------------------------------------------------
#                     ANNOTATIONS FOR REFERENCE AND QUERY
##---------------------------------------------------------------------------

# Read annotation data
ref.annotation <- read.delim("annotation.txt", sep = "\t", check.names = FALSE)
ref.annotation = ref.annotation %>% remove_rownames %>% 
  column_to_rownames(var = "Sample")

# Read table design data and modify it
query.annotation <- read.delim("annotation_paper.txt", sep = "\t", 
                          check.names = FALSE)
query.annotation = query.annotation[, -3]
query.annotation = query.annotation[, -3]
query.annotation = query.annotation[-29,]
query.annotation = query.annotation %>% remove_rownames %>% 
  column_to_rownames(var = "sample")
query.annotation=as.data.frame(query.annotation)

##---------------------------------------------------------------------------

##================================================================
##                   SECTION 2: Integration                     ==
##================================================================


# Create Seurat objects for reference data
reference = CreateSeuratObject(
  counts = local.cell.line,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = ref.annotation
)

# Create Seurat objects for query data
query = CreateSeuratObject(
  counts = data,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = query.annotation,
  project = "my_project2"
)
query <- SetIdent(query, value = query@meta.data$group)

# Pre-process the reference dataset without integration
reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference, approx = FALSE)
reference <- FindNeighbors(reference, dims = 1:4)
reference <- FindClusters(reference)

reference <- RunUMAP(reference, dims = 1:4,reduction = "pca")

DimPlot(reference, group.by = c("Cell_type"), label = TRUE)


# Normalize peptide data
query <- NormalizeData(query)

# Find anchors for integration
anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:6, 
                               reference.reduction = "pca",k.anchor=7,
                               k.score = 27)
predictions <- TransferData(anchorset = anchors, refdata = reference$Cell_type, 
                            dims = 1:4, k.weight = 8)

# Add predicted cell types to the peptide dataset
ovarian.query <- AddMetaData(query, metadata = predictions)

# Evaluate the quality of the predictions
ovarian.query$prediction.match <- 
  ovarian.query$predicted.id == ovarian.query$group
table(ovarian.query$prediction.match)

# Display the distribution of predicted cell types
table(ovarian.query$predicted.id)

# Run UMAP on the reference dataset for visualization
ovarian.ref <- RunUMAP(reference, dims = 1:6, return.model = TRUE)

# Map the query dataset to the reference UMAP
ovarian.query <- MapQuery(anchorset = anchors, reference = ovarian.ref, 
                          query = ovarian.query, 
                          refdata = list(predicted_celltype = "Cell_type"),
                          reference.reduction = "pca", reduction.model = "umap",
                          transferdata.args = list(k.weight = 8))



# Plot side-by-side UMAPs of the reference and query datasets
p1 <- DimPlot(ovarian.ref, reduction = "umap", group.by = "Cell_type", 
              label = TRUE, label.size = 5, pt.size = 5, repel = TRUE) + 
  ggtitle("Reference annotations")
p2 <- DimPlot(ovarian.query, reduction = "ref.umap", 
              group.by = "predicted.predicted_celltype", label = TRUE,
              label.size = 5, pt.size = 5, repel = TRUE, alpha=1) + 
  ggtitle("Query transferred labels")
p1+p2
p3 <- DimPlot(ovarian.query,pt.size =4,label=TRUE,repel = TRUE)+
  ggtitle("Query annotations")
plot=p1+p3+p2+plot_annotation(title = "Proteome_VS_Proteome")  &  
  theme(plot.title = element_text(hjust = 0.5,size = 15)) 
plot
