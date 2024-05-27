# Load necessary libraries
library(Seurat)
library(SeuratObject)
#library(SeuratData)
library(ggplot2)
library("writexl")
library(tidyverse)
library(patchwork)

# Load the data into R
data <- read.delim("Phospho.pg_matrix(protein level).tsv", sep = "\t", check.names = FALSE)

# Remove unnecessary columns from the data
data = subset(data, select = -c(3:5))
data = subset(data, select = -c(1))

# Extract gene names and transpose the data
Protein = data[, 1]
rownames(data) = make.names(Protein, unique = TRUE)
data = data[, -1]
data = as.data.frame(data)
#datat = apply(data, 2, median, na.rm = TRUE)
#data=2^data
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

# Read peptide data and preprocess it
peptide <- read.delim("abundance_protein_MD.tsv", sep = "\t", check.names = FALSE)
peptide <- peptide[, !(names(peptide) %in% c("NumberPSM", "Gene", "MaxPepProb", "ReferenceIntensity"))]
Gene = peptide[, 1]

# Make unique Gene names
rownames(peptide) = make.names(Gene, unique = TRUE)

# Delete column 1
peptide = peptide[, -1]
peptide = 2^peptide

#select only complete cases
#peptide <- peptide[complete.cases(peptide), ]


# Let's say 'df' is your data frame and you want to remove columns 'col1' and 'col2'
#cols_to_remove <- c("JHUQC", "JHUQC.1","JHUQC.2","JHUQC.3","JHUQC.4")

# Remove the columns
#peptide <- peptide[, !(names(peptide) %in% cols_to_remove)]

##################### for each SAMPLE############################
# Impute low value for missing data in peptide
# Identify numeric columns
numeric_cols <- sapply(peptide, is.numeric)

## Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(peptide[, i], na.rm = TRUE)
  
  #  # Define a low value (e.g., min_val - 1)
  low_val <- min_val
  
  #  # Replace NA values with the low value
  peptide[is.na(peptide[, i]), i] <- low_val
}

##################### ##########################################


###########################FOR EACH GENE#########################
# Assuming data is your data frame

# Identify numeric columns
#numeric_cols <- sapply(peptide, is.numeric)

# Iterate over rows
#for (row_index in 1:nrow(peptide)) {
# Extract values for the current row
#  row_values <- peptide[row_index, numeric_cols]

# Calculate the minimum value for the row, excluding NA values
#  min_val <- min(row_values, na.rm = TRUE)

# Replace missing values in the current row with the calculated minimum
#  peptide[row_index, numeric_cols][is.na(row_values)] <- min_val
#}

##################### ##########################################

# Read annotation data
annotation <- read.delim("annotation.txt", sep = "\t", check.names = FALSE)
ann = annotation %>% remove_rownames %>% column_to_rownames(var = "Sample")

# Read table design data and modify it
annotation2 <- read.delim("table_design.tsv", sep = "\t", check.names = FALSE)
annotation2$group <- gsub("_", "", annotation2$group)
annotation2$new_sample <- paste(annotation2$group, annotation2$sample, sep = "_")

# Create a mapping between old and new column names
column_mapping <- setNames(annotation2$new_sample, annotation2$sample)

# Rename the columns in the peptide dataset
colnames(peptide) <- column_mapping[colnames(peptide)]

annotation2 = annotation2[, -1]
ann2 = annotation2 %>% remove_rownames %>% column_to_rownames(var = "new_sample")

# Create Seurat objects for RNA data
sdata = CreateSeuratObject(
  counts = data,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = ann,
  project = "my_project"
)

# Create Seurat objects for peptide data
peptides = CreateSeuratObject(
  counts = peptide,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = ann2,
  project = "my_project2"
)

# Pre-process the RNA dataset without integration
sdata <- NormalizeData(sdata)
sdata <- FindVariableFeatures(sdata)
sdata <- ScaleData(sdata)
sdata <- RunPCA(sdata, approx = FALSE)
sdata <- FindNeighbors(sdata, dims = 1:6)
sdata <- FindClusters(sdata)
head(Idents(sdata), 54)
sdata <- RunUMAP(sdata, dims = 1:5,reduction = "pca")
DimPlot(sdata, group.by = c("Cell_type"), label = TRUE)

# Normalize peptide data
peptides <- NormalizeData(peptides)

# Find anchors for integration
anchors <- FindTransferAnchors(reference = sdata, query = peptides, dims = 1:4, reference.reduction = "pca",k.anchor=10)
predictions <- TransferData(anchorset = anchors, refdata = sdata$Cell_type, dims = 1:4, k.weight = 6)

# Add predicted cell types to the peptide dataset
pancreas.query <- AddMetaData(peptides, metadata = predictions)

# Evaluate the quality of the predictions
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$group
table(pancreas.query$prediction.match)

# Display the distribution of predicted cell types
table(pancreas.query$predicted.id)

# Visualize expression of specific genes in predicted cell types
#VlnPlot(pancreas.query, c("TP53", "BRCA1", "BRCA2", "PIK3CA", "PTEN", "MYC", "ARID1A", "ERBB2", "BRAF", "PALB2", "CHEK2", "BARD1", "MRE11", "CTNNB1", "PTEN"), group.by = "predicted.id")

# Run UMAP on the reference dataset for visualization
pancreas.ref <- RunUMAP(sdata, dims = 1:4, return.model = TRUE)

# Map the query dataset to the reference UMAP
pancreas.query <- MapQuery(anchorset = anchors, reference = pancreas.ref, query = pancreas.query, refdata = list(predicted_celltype = "Cell_type"),
                           reference.reduction = "pca", reduction.model = "umap", transferdata.args = list(k.weight =6))

# Plot side-by-side UMAPs of the reference and query datasets
p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "Cell_type", label = TRUE, label.size = 5, pt.size = 5, repel = TRUE) + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.predicted_celltype", label = TRUE,
              label.size = 5, pt.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
#p1 + p2
p3 <- DimPlot(pancreas.query,pt.size = 1)+ggtitle("Query annotations")
plot=p1+p3+p2+plot_annotation(title = "Phospo_VS_Phospo on Protein Level") &  theme(plot.title = element_text(hjust = 0.5,size = 15)) 
plot
#select.cells <- CellSelector(plot = p3)
#select.cells

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
# #merge reference and query
# reference$id <- 'reference'
# pbmc3k$id <- 'query'
# refquery <- merge(reference, pbmc3k)
# refquery[["umap"]] <- merge(reference[["umap"]], pbmc3k[["ref.umap"]])
# refquery <- RunUMAP(refquery, reduction = 'umap', dims = 1:4)
# p4=DimPlot(refquery, group.by = 'id')+ggtitle("Reference and Query annotations")
# plot2=p1+p3+p4+p2+plot_annotation(title = "Phospo_VS_Phospo on on Protein Level")  &  theme(plot.title = element_text(hjust = 0.5,size = 15)) 
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


# VlnPlot(pancreas.query, features = c("MYH11-S1954", "GPBP1L1-S445"),group.by = "predicted.id")
# FeaturePlot(pancreas.query, features = c("MYH11-S1954"))
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



all.markers <- FindAllMarkers(object = pancreas.query)
all.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),group_by=TRUE) %>% 
  top_n(avg_log2FC,n=10) %>% 
  ungroup() -> top10
pancreas.query <- ScaleData(pancreas.query)
pancreas.query$Cells=Cells(pancreas.query)
P1=DoHeatmap(pancreas.query,features = top10$gene,group.bar = TRUE,label=TRUE,angle=30,size = 3.5)+  guides(colour=FALSE)
P1
# P2=DoHeatmap(pancreas.query, features = top10$gene,group.by = "Cells",group.bar = TRUE,label=TRUE,size = 3,angle=90)+  guides(colour=FALSE)
# P2
# P1+P2




all.markers2 <- FindAllMarkers(object = pancreas.ref)
head(x = all.markers2)
top10 <- all.markers2 %>% group_by(cluster) %>% top_n(20, avg_log2FC)
all.markers2 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(pancreas.ref, features = top10$gene,group.bar = TRUE)


names(all.markers) <- paste0(names(all.markers), ".query")
all.markers$gene <- rownames(all.markers)
names(all.markers2) <- paste0(names(all.markers2), ".reference")
all.markers2$gene <- rownames(all.markers2)

merge_dat <- merge(all.markers2, all.markers, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.query), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.query < 0.005 & 
                                 merge_dat$p_val.reference < 0.005)]
only_reference <- merge_dat$gene[which(merge_dat$p_val.query > 0.005 & 
                                         merge_dat$p_val.reference < 0.005)]
only_query <- merge_dat$gene[which(merge_dat$p_val.query < 0.005 & 
                                     merge_dat$p_val.reference > 0.005)]
print(paste0('# Common: ',length(common)))
print(paste0('# Only in reference: ',length(only_reference)))
print(paste0('# Only in query: ',length(only_query)))




# pbmc <- FindVariableFeatures(pancreas.query, selection.method = "vst", nfeatures = 46275)
# top10 <- head(VariableFeatures(pbmc), 10)
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2



# add a column of NAs
all.markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
all.markers$diffexpressed[all.markers$avg_log2FC.query > 1 & all.markers$p_val.query < 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
all.markers$diffexpressed[all.markers$avg_log2FC.query < -1 & all.markers$p_val.query < 0.5] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=all.markers, aes(x=avg_log2FC.query, y=-log10(p_val.query), col=diffexpressed)) + geom_point() + theme_minimal()
p
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")
p2

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
p3
# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
all.markers$delabel <- NA
all.markers$delabel[all.markers$diffexpressed != "NO"] <- all.markers$gene.query[all.markers$diffexpressed != "NO"]
library(ggrepel)
ggplot(data=all.markers, aes(x=avg_log2FC.query, y=-log10(p_val.query), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel( max.overlaps = getOption("ggrepel.max.overlaps", default = 3)) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")





all.markers$delabel <- NA
all.markers$delabel[all.markers$diffexpressed != "NO"] <- all.markers$cluster.query[all.markers$diffexpressed != "NO"]
ggplot(data=all.markers, aes(x=avg_log2FC.query, y=-log10(p_val.query), col=diffexpressed, label=gene.query)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()
