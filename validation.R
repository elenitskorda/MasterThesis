# Load necessary libraries
library(Seurat)
library(SeuratObject)
#library(SeuratData)
library(ggplot2)
library("writexl")
library(tidyverse)
library(patchwork)
library(enrichR)
# Load the data into R
data <- read.delim("reference_paper.txt", sep = "\t", check.names = FALSE)
#"OvCa Proteins 5sub Loess DE_stats.tsv"


# Extract gene names and transpose the data
Gene = data[, 1]
rownames(data) = make.names(Gene, unique = TRUE)   
data = data[, -1]
data = as.data.frame(data)
#datat = apply(data, 2, median, na.rm = TRUE)
data=2^data
#data <- data[complete.cases(data), ]
# Impute median value for missing data

##################### for each SAMPLE############################
# Identify numeric columns
numeric_cols <- sapply(data, is.numeric)


# Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(data[, i], na.rm = TRUE)
  
  # Define a low value (e.g., min_val - 1)
  low_val <- min_val
  
  #  # Replace NA values with the low value
  data[is.na(data[, i]), i] <- low_val
}
##################### ##########################################



###########################FOR EACH GENE#########################
# Assuming data is your data frame

# Identify numeric columns
#numeric_cols <- sapply(data, is.numeric)

# Iterate over rows
#for (row_index in 1:nrow(data)) {
# Extract values for the current row
#  row_values <- data[row_index, numeric_cols]

# Calculate the minimum value for the row, excluding NA values
#  min_val <- min(row_values, na.rm = TRUE)

# Replace missing values in the current row with the calculated minimum
#  data[row_index, numeric_cols][is.na(row_values)] <- min_val
#}

##################### ##########################################

# Read peptide peptide and preprocess it
peptide <- read.delim("OvCa Proteins 5sub Loess DE_stats.tsv", sep = "\t", check.names = FALSE)
# Remove unnecessary columns from the peptide
peptide = subset(peptide, select = -c(1:3))
peptide = subset(peptide, select = -c(2:12))

# Extract gene names and transpose the peptide
Gene = peptide[, 1]
rownames(peptide) = make.names(Gene, unique = TRUE)
peptide = peptide[, -1]

#peptidet = apply(peptide, 2, median, na.rm = TRUE)
peptide=2^peptide
#peptide <- peptide[complete.cases(peptide), ]
# Impute median value for missing peptide

##################### for each SAMPLE############################
# Identify numeric columns
numeric_cols <- sapply(peptide, is.numeric)


# Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(peptide[, i], na.rm = TRUE)
  
  # Define a low value (e.g., min_val - 1)
  low_val <- min_val
  
  #  # Replace NA values with the low value
  peptide[is.na(peptide[, i]), i] <- low_val
}




##################### ##########################################

# Read annotation data
annotation <- read.delim("annotation.txt", sep = "\t", check.names = FALSE)
ann = annotation %>% remove_rownames %>% column_to_rownames(var = "Sample")

# Read table design data and modify it
annotation2 <- read.delim("annotation_paper.txt", sep = "\t", check.names = FALSE)
annotation2 = annotation2[, -3]
annotation2 = annotation2[, -3]
annotation2 = annotation2[-29,]
ann2 = annotation2 %>% remove_rownames %>% column_to_rownames(var = "sample")
ann2=as.data.frame(ann2)

# Create Seurat objects for RNA data
sdata = CreateSeuratObject(
  counts = peptide,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = ann
)

# Create Seurat objects for peptide data
peptides = CreateSeuratObject(
  counts = data,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = ann2,
  project = "my_project2"
)
peptides <- SetIdent(peptides, value = peptides@meta.data$group)
# Pre-process the RNA dataset without integration
sdata <- NormalizeData(sdata)
sdata <- FindVariableFeatures(sdata)
sdata <- ScaleData(sdata)
sdata <- RunPCA(sdata, approx = FALSE)
sdata <- FindNeighbors(sdata, dims = 1:4)
sdata <- FindClusters(sdata)

sdata <- RunUMAP(sdata, dims = 1:4,reduction = "pca")

DimPlot(sdata, group.by = c("Cell_type"), label = TRUE)


# Normalize peptide data
peptides <- NormalizeData(peptides)

# Find anchors for integration
anchors <- FindTransferAnchors(reference = sdata, query = peptides, dims = 1:6, reference.reduction = "pca",k.anchor=7,k.score = 27)
predictions <- TransferData(anchorset = anchors, refdata = sdata$Cell_type, dims = 1:4, k.weight = 8)

# Add predicted cell types to the peptide dataset
pancreas.query <- AddMetaData(peptides, metadata = predictions)

# Evaluate the quality of the predictions
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$group
table(pancreas.query$prediction.match)

# Display the distribution of predicted cell types
table(pancreas.query$predicted.id)


# Run UMAP on the reference dataset for visualization
pancreas.ref <- RunUMAP(sdata, dims = 1:6, return.model = TRUE)

# Map the query dataset to the reference UMAP
pancreas.query <- MapQuery(anchorset = anchors, reference = pancreas.ref, query = pancreas.query, refdata = list(predicted_celltype = "Cell_type"),
                           reference.reduction = "pca", reduction.model = "umap", transferdata.args = list(k.weight = 8))



# Plot side-by-side UMAPs of the reference and query datasets
p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "Cell_type", label = TRUE, label.size = 5, pt.size = 5, repel = TRUE) + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.predicted_celltype", label = TRUE,
              label.size = 5, pt.size = 3, repel = TRUE, alpha=2) + ggtitle("Query transferred labels")
p1+p2
p3 <- DimPlot(pancreas.query,pt.size =3,label=TRUE)+ggtitle("Query annotations")
plot=p1+p3+p2+plot_annotation(title = "Proteome_VS_Proteome")  &  theme(plot.title = element_text(hjust = 0.5,size = 15)) 
plot
#select.cells <- CellSelector(plot = p3)
#head(select.cells)

# 
# tumor_markers <- FindMarkers(pancreas.query, ident.1 = "Tumor", only.pos = TRUE, logfc.threshold = 0.1)
# print(head(tumor_markers))
# 
# normal_markers <- FindMarkers(pancreas.query, ident.1 = "Normal", only.pos = TRUE, logfc.threshold = 0.1)
# print(head(normal_markers))
# 
# 
# reference <- DietSeurat(pancreas.ref, counts = FALSE, dimreducs = "umap")
# 
# pbmc3k <- DietSeurat(pancreas.query, counts = FALSE, dimreducs = "ref.umap")
# 
# 
# # #merge reference and query
# reference$id <- 'reference'
# pbmc3k$id <- 'query'
# refquery <- merge(reference, pbmc3k)
# refquery[["umap"]] <- merge(reference[["umap"]], pbmc3k[["ref.umap"]])
# refquery <- RunUMAP(refquery, reduction = 'umap', dims = 1:4)
# p4=DimPlot(refquery, group.by = 'id')+ggtitle("Reference and Query annotations")
# plot2=p1+p3+p4+p2+plot_annotation(title = "Proteome_VS_Proteome")  &  theme(plot.title = element_text(hjust = 0.5,size = 15)) 
# plot2
# 
# pbmc <- FindVariableFeatures(pancreas.query, selection.method = "vst", nfeatures = 46275)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# 
# 
# # VlnPlot(pancreas.query, features = c("MYH11-S1954", "GPBP1L1-S445"),group.by = "predicted.id")
# # FeaturePlot(pancreas.query, features = c("MYH11-S1954"))
# 
# pbmc.markers <- FindAllMarkers(pancreas.query, only.pos = TRUE)
# pbmc.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1)
# 
# pbmc.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 10) %>%
#   ungroup() -> top10
# pancreas.query <- ScaleData(pancreas.query)
# DoHeatmap(pancreas.query, features = top10$gene)

###########################################add from phospho vs phospho#########################################################



Idents(pancreas.query)<-pancreas.query$predicted.id
Idents(pancreas.query)

#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')
library(presto)




#merge Reference and query

reference = DietSeurat(pancreas.ref, counts = FALSE, data = TRUE,scale.data = FALSE,dimreducs = "umap")
Idents(reference)=reference$Cell_type

# Idents(reference, cells = 1:6) <- 'RefA2780'
# Idents(reference, cells = 7:12) <- 'RefCAOV3'
# Idents(reference, cells = 13:18) <- 'RefES2'
# Idents(reference, cells = 19:24) <- 'RefIGROV1'
# Idents(reference, cells = 25:30) <- 'RefJHOC5'
# Idents(reference, cells = 31:36) <- 'RefKURAM'
# Idents(reference, cells = 37:42) <- 'RefNIHOVCAR3'
# Idents(reference, cells = 43:48) <- 'RefOVSAHO'
# Idents(reference, cells = 49:54) <- 'RefSKOV3' 
#ScaleData(pancreas.query)
pbmc3k = DietSeurat(pancreas.query, counts = FALSE, data = TRUE,scale.data = FALSE,dimreducs = "ref.umap",assays="RNA")
#pbmc3k <- pancreas.query


reference$id <- 'reference'
pbmc3k$id <- 'query'
refquery <- merge(reference, y=pbmc3k)



refquery[["umap"]] <- merge(reference[["umap"]], pbmc3k[["ref.umap"]])

# install.packages('BiocManager')
# BiocManager::install('multtest')
# install.packages('metap')

library(multtest)
library(metap)
DefaultAssay(refquery) <- "RNA"
refquery[['groups']] <-refquery$id

refquery=JoinLayers(refquery)
#refquery=DietSeurat(refquery, counts = TRUE, data = TRUE,scale.data = FALSE,dimreducs = "umap")

conserved=FindConservedMarkers(object =refquery,ident.1 = 'KURAM',grouping.var = "groups")



all.markers <- FindAllMarkers(object = pancreas.query,only.pos = FALSE)

write.csv(all.markers, file="all_markers_proteome.csv")
all.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),group_by=TRUE) %>% 
  top_n(avg_log2FC,n=10) %>% 
  ungroup() -> top10

pancreas.query <- ScaleData(pancreas.query)
pancreas.query$Cells=Cells(pancreas.query)
P1=DoHeatmap(pancreas.query,features = top10$gene,group.bar = TRUE,label=TRUE,angle=30,size = 3.5)+  guides(colour=FALSE)
P1
# P2=DoHeatmap(pancreas.query, features = top10$gene,group.by = "Cells",group.bar = TRUE,label=TRUE,size = 2.5,angle=90)+  guides(colour=FALSE)
# P2
# P1+P2

#Import table with common genes between this and pathway analysis
common_genes<-read.csv("pathway_exp_commongenes.csv")
P2=DoHeatmap(pancreas.query,features = common_genes$Genes,group.bar = TRUE,label=TRUE,angle=30,size = 3.5)+  guides(colour=FALSE)
P2

all.markers %>%
  filter(gene=="MTDH" | gene=="SND1")

