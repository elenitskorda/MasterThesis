# -*- coding: utf-8 -*-

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Date: Mon Sep 5 10:28:18 2023                                                                        ::
##  File Name: phosphoproteome_integration_protein.R                                                     ::
##  Author: Eleni Theofania Skorda                                                                       ::
##                                                                                                       ::
##  Description:                                                                                         ::
##      Integrates phospoproteomic local cell line data (Reference) with CPTAC data (Query)              ::
##      in order to indetify cell lines within the query dataset.                                        ::
##                                                                                                       ::
##  List of Functions:                                                                                   ::
##      none                                                                                             ::
##                                                                                                       ::
##  Procedure:                                                                                           ::
##      1. Prepare reference and query datasets to a proper format                                       ::
##      2. Performing integration                                                                        ::
##      3. Finds markers (differentially expressed genes) for each of the identity classes in a dataset  ::
##                                                                                                       ::
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library("writexl")
library(tidyverse)
library(patchwork)
library(presto)
library(multtest)
library(metap)

##==========================================================================
##  SECTION 1: Prepare reference and query datasets to a proper format    ==
##==========================================================================


#                            REFERENCE DATASET
##---------------------------------------------------------------------------

# Load the local.cell.lines into R
local.cell.lines <- read.delim("Phospho.pg_matrix(protein level).tsv", 
                               sep = "\t", check.names = FALSE)

# Remove unnecessary columns from the local.cell.lines
local.cell.lines = subset(local.cell.lines, select = -c(3:5))
local.cell.lines = subset(local.cell.lines, select = -c(1))

# Extract gene names and transpose the local.cell.lines
Protein = local.cell.lines[, 1]
rownames(local.cell.lines) = make.names(Protein, unique = TRUE)
local.cell.lines = local.cell.lines[, -1]
local.cell.lines = as.data.frame(local.cell.lines)

# Impute median value for missing data
# Identify numeric columns
numeric_cols <- sapply(local.cell.lines, is.numeric)


# Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(local.cell.lines[, i], na.rm = TRUE)
  
  # Define a low value (e.g., min_val - 1)
  low_val <- min_val
  
  #  # Replace NA values with the low value
  local.cell.lines[is.na(local.cell.lines[, i]), i] <- low_val
}

##---------------------------------------------------------------------------
#                              QUERY DATASET
##---------------------------------------------------------------------------

# Read CPTAC data and preprocess it
CPTAC <- read.delim("abundance_protein_MD.tsv", sep = "\t", check.names = FALSE)
CPTAC <- CPTAC[, !(names(CPTAC) %in% c("NumberPSM", 
                                       "Gene", 
                                       "MaxPepProb", 
                                       "ReferenceIntensity"))]
Gene = CPTAC[, 1]

# Make unique Gene names
rownames(CPTAC) = make.names(Gene, unique = TRUE)

# Delete column 1
CPTAC = CPTAC[, -1]
CPTAC = 2^CPTAC

# Impute low value for missing data in CPTAC
# Identify numeric columns
numeric_cols <- sapply(CPTAC, is.numeric)

## Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(CPTAC[, i], na.rm = TRUE)
  
  #  # Define a low value 
  low_val <- min_val
  
  #  # Replace NA values with the low value
  CPTAC[is.na(CPTAC[, i]), i] <- low_val
}

##---------------------------------------------------------------------------
#                     ANNOTATIONS FOR REFERENCE AND QUERY
##---------------------------------------------------------------------------

# Read annotation data
ref.annotation <- read.delim("annotation.txt", sep = "\t", check.names = FALSE)
ref.annotation = ref.annotation %>% remove_rownames %>% 
  column_to_rownames(var = "Sample")

# Read table design data and modify it
query.annotation <- read.delim("table_design.tsv", sep = "\t", 
                               check.names = FALSE)
query.annotation$group <- gsub("_", "", query.annotation$group)
query.annotation$new_sample <- paste(query.annotation$group, 
                                     query.annotation$sample, sep = "_")

# Create a mapping between old and new column names
column_mapping <- setNames(query.annotation$new_sample, query.annotation$sample)

# Rename the columns in the query dataset
colnames(CPTAC) <- column_mapping[colnames(CPTAC)]

query.annotation = query.annotation[, -1]
query.annotation = query.annotation %>% remove_rownames %>% 
  column_to_rownames(var = "new_sample")

##---------------------------------------------------------------------------

##================================================================
##                   SECTION 2: Integration                     ==
##================================================================

# Create Seurat objects for RNA data
reference = CreateSeuratObject(
  counts = local.cell.lines,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = ref.annotation,
  project = "my_project"
)

# Create Seurat objects for peptide data
query = CreateSeuratObject(
  counts = CPTAC,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = query.annotation,
  project = "my_project2"
)

# Pre-process the RNA dataset without integration
reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference, approx = FALSE)
reference <- FindNeighbors(reference, dims = 1:6)
reference <- FindClusters(reference)
head(Idents(reference), 54)
reference <- RunUMAP(reference, dims = 1:5,reduction = "pca")
DimPlot(reference, group.by = c("Cell_type"), label = TRUE)

# Normalize peptide data
query <- NormalizeData(query)

# Find anchors for integration
anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:4, 
                               reference.reduction = "pca",k.anchor=10)
predictions <- TransferData(anchorset = anchors, refdata = reference$Cell_type,
                            dims = 1:4, k.weight = 6)

# Add predicted cell types to the peptide dataset
ovarian.query <- AddMetaData(query, metadata = predictions)

# Evaluate the quality of the predictions
ovarian.query$prediction.match <- ovarian.query$predicted.id == 
  ovarian.query$group
table(ovarian.query$prediction.match)

# Display the distribution of predicted cell types
table(ovarian.query$predicted.id)

# Run UMAP on the reference dataset for visualization
ovarian.ref <- RunUMAP(reference, dims = 1:4, return.model = TRUE)

# Map the query dataset to the reference UMAP
ovarian.query <- MapQuery(anchorset = anchors, reference = ovarian.ref, 
                          query = ovarian.query, 
                          refdata = list(predicted_celltype = "Cell_type"),
                          reference.reduction = "pca", reduction.model = "umap",
                          transferdata.args = list(k.weight =6))

# Plot side-by-side UMAPs of the reference and query datasets
p1 <- DimPlot(ovarian.ref, reduction = "umap", group.by = "Cell_type", 
              label = TRUE, label.size = 5, pt.size = 5, repel = TRUE) + 
  ggtitle("Reference annotations")
p2 <- DimPlot(ovarian.query, reduction = "ref.umap", 
              group.by = "predicted.predicted_celltype", label = TRUE,
              label.size = 5, pt.size = 3, repel = TRUE) + 
  ggtitle("Query transferred labels")
#p1 + p2
p3 <- DimPlot(ovarian.query,pt.size = 1)+ggtitle("Query annotations")
plot=p1+p3+p2+plot_annotation(title = "Phospo_VS_Phospo on Protein Level") &
  theme(plot.title = element_text(hjust = 0.5,size = 15)) 
plot


##================================================================
##  SECTION 3: Finds markers (differentially expressed genes)   ==
##            for each of the identity classes in a dataset     ==
##================================================================

# Rename the idents of ovarian.query with integrated ones
Idents(ovarian.query)<-ovarian.query$predicted.id
Idents(ovarian.query)


# Merge Reference and query
# Keep only certain aspects of the Seurat object. Useful for merge()
reference = DietSeurat(ovarian.ref, counts = FALSE, data = TRUE,
                       scale.data = FALSE,dimreducs = "umap")
Idents(reference)=reference$Cell_type
query = DietSeurat(ovarian.query, counts = FALSE, data = TRUE,
                   scale.data = FALSE,dimreducs = "ref.umap",assays="RNA")
reference$id <- 'reference'
query$id <- 'query'

# Merge the reference and query Seurat Obj
reference.query <- merge(reference, y=query)


# Merge the reductions
reference.query[["umap"]] <- merge(reference[["umap"]], query[["ref.umap"]])
DefaultAssay(reference.query) <- "RNA"
reference.query[['groups']] <-reference.query$id

# Split and Join Layers Together
reference.query=JoinLayers(reference.query)


# Finds differentially expressed genes for each of the identity classes in a dataset
all.markers <- FindAllMarkers(object = ovarian.query)

# Group the markers by cluster & select top 10 markers for each cluster
all.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),group_by=TRUE) %>% 
  top_n(avg_log2FC,n=10) %>% 
  ungroup() -> top10

# Scales and centers features in the query dataset
ovarian.query <- ScaleData(ovarian.query)
ovarian.query$Cells=Cells(ovarian.query)

# Plot top 10 markers for each cluster in heatmap
P1=DoHeatmap(ovarian.query,features = top10$gene,group.bar = TRUE,
             label=TRUE,angle=30,size = 3.5)+  guides(colour=FALSE)

# After running the above "DoHeatmap" it might appear the following warning

##-------------------------------------------------------------------------------------------
##  Warning message: The `<scale>` argument of `guides()` cannot be `FALSE`.               --
##  Use "none" instead as of ggplot2 3.3.4. This warning is displayed once every 8 hours.  --
##  Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.   --
##-------------------------------------------------------------------------------------------

# IF YES RUN IT AGAIN
P1

