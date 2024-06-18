# -*- coding: utf-8 -*-



# Load necessary libraries
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(factoextra)
library(limma)
library(ggpubr)
library(PCAtools)
library(ggplot2)

# Import ratio peptide data from CPTAC
CPTAC <- read.delim("ratio_peptide_MD.tsv", sep = "\t",check.names = FALSE)

# Import plex annotation from CPTAC
plex <- read.delim("experiment_annotation2.txt", sep = "\t",check.names = FALSE)

# Import clinical annotation file from CPTAC
clinical = read.delim("PDC_clinical_manifest_09042023_135706.tsv", sep = "\t",
                      check.names = FALSE)

# Import samples annotation file from CPTAC
samples <- read.delim("table_design.tsv", sep = "\t",check.names = FALSE)

# Import clinical annotation file from CPTAC
clinical2 = read.delim("clinical2.txt", sep = "\t",check.names = FALSE)

# Curate the dataset by subsetting columns 1 and 4
clinical2 = subset(clinical2, select = c(1, 4))
sample=clinical2[,1]

# Make samples name unique
rownames(clinical2) = make.names(sample, unique=TRUE)

# Replace X with " "
new <- gsub("X", "", rownames(clinical2))

rownames(clinical2)=new
clinical2=clinical2[, -1]

# Delete uneccessary columns from CPTAC data
CPTAC <- CPTAC[, !(names(CPTAC) %in% c("Index", 
                                       "ProteinID",
                                       "Peptide",
                                       "MaxPepProb",
                                       "ReferenceIntensity"))]

# Remove NAs
CPTAC=na.omit(CPTAC)

# Name first column as "Gene"
Gene=CPTAC[,1]

#Make unique Gene names
rownames(CPTAC) = make.names(Gene, unique=TRUE)

#Delete column 1
CPTAC=CPTAC[,-1]


#Calculate the variance for each gene in a column called rowVar
CPTAC$row_var = rowVars(as.matrix(CPTAC[,c(-1)]))

# Check the range of row var
range <- range(CPTAC$row_var)
row_var=CPTAC$row_var

#Sort the values in a decreasing order
sorted_values <- sort(row_var, decreasing = TRUE)

#set an index of 50%
split_index <- length(sorted_values) * 0.5

#Choose the top 50 %
top_50_percent <- sorted_values[1:split_index]
keep_rows <- row_var %in% top_50_percent
filtered_df <- CPTAC[keep_rows, ]

# Transpose the filtered data frame
CPTAC=t(filtered_df)

# Name column 1 as sample
sample=samples[,1]

# Subset columns from clinical file 2,8,12,14
clinical = subset(clinical, select = c(2, 8, 12, 14))

# Define column names in clinical annotation
colnames(clinical)[1] ="sample"
colnames(clinical) <- c("sample","Race",'TissueOrgan',"TumorStage")

# Define sample names in clinical annotation
clinical$sample[clinical$sample == '01OV007'] <- '01OV007.1'
clinical$sample[clinical$sample == '01OV017'] <- '01OV017.1'
clinical$sample[clinical$sample == '01OV023'] <- '01OV023.1'
clinical$sample[clinical$sample == '17OV001'] <- '17OV001.1'
clinical$sample[clinical$sample == '01OV039'] <- '01OV039.1'
clinical$sample[clinical$sample == '11OV002'] <- '11OV002.1'
clinical$sample[clinical$sample == '14OV011'] <- '01OV011.1'

 new_row <- data.frame(
   sample = c('01OV029.1','15OV001.1','01OV007','01OV017',
              '01OV023','17OV001','01OV039','11OV002','17OV002.1','14OV011'),
   Race = c("No Data","No Data","No Data","No Data","No Data",
            "No Data","No Data","No Data","No Data","No Data"),
   TissueOrgan= c("No Data","No Data","No Data","No Data","No Data",
                  "No Data","No Data","No Data","No Data","No Data"),
   TumorStage=c("No Data","No Data","No Data","No Data","No Data",
                "No Data","No Data","No Data","No Data","No Data"))
 clinical <- rbind(clinical, new_row)

# Merge samples and clinical annotations by sample
predesign=merge(samples, clinical,by="sample")

# Merge predesign and plex annotations by sample
design=merge(predesign, plex,by="sample")

# Subset columns from design annotations
design=subset(design,select=c(1:6))

# Name the subseted columns
colnames(design) <- c("sample",'HealthStatus',"Race",'Organ',
                      "TumorStage", "Plex")

# Make samples unique samples
row.names(design) <- make.names(design$sample, unique=TRUE)

# Replace "x" from samples names with " "
new <- gsub("X", "", rownames(design))

rownames(design)=new
design=design[,-1]

# Check id rownames matches between annotation and dataset
row_names_match <- all(rownames(design) == rownames(CPTAC))

# Check which names are not matchinf
row_names_df1 <- rownames(CPTAC)
row_names_df2 <- rownames(design)
non_matching_row_names <- setdiff(row_names_df1, row_names_df2)

# Preserve the common rows between annotation file and dataset
common_row_names <- intersect(row_names_df1, row_names_df2)

# Filter both data frames to keep only the common row names
CPTAC <- CPTAC[row_names_df1 %in% common_row_names, ]

# Double check again if the row names are matching
row_names_df1 <- rownames(CPTAC)
row_names_df2 <- rownames(design)
non_matching_row_names <- setdiff(row_names_df1, row_names_df2)


if (length(non_matching_row_names) == 0) {
  # All row names match between df1 and df2
  cat("All row names match between df1 and df2.\n")
} else {
  cat("Row names that do not match between df1 and df2:\n")
  print(non_matching_row_names)
}
# All row names match between df1 and df2

# Make the annotation file as matrix
as.matrix(design)

# Reorder rows of CPTAC to match the annotation file
CPTAC <- CPTAC[rownames(design), ]

# Transpose CPTAC dataset
CPTAC=t(CPTAC)

# Make row names of annotation file same as the colnames from CPTAC dataset
row.names(design) <- colnames(CPTAC)

#PCATOOLS


# Implement Principal component Analysis
p <- pca(CPTAC,metadata = design)

# Correlate principal components to continuous variable metadata and test significancies of these.
eigencorplot(p,components = getComponents(p, 1:10),
             metavars = c('HealthStatus','Race','Organ','TumorStage','Plex'))

# Perform Horn's parallel analysis to choose the number of principal components to retain.
horn <- parallelPCA(CPTAC)
horn$n

# Find the elbow point in the curve of variance explained by each successive PC.
# This can be used to determine the number of PCs to retain.
elbow <- findElbowPoint(p$variance)
elbow


screeplot(p,
          components = getComponents(p, 1:30),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

# Generate plots

# Showing the distribution of explained variance across all or select principal components 
pscree <- screeplot(p, components = getComponents(p, 1:35),
                    hline = 80, vline = 33, axisLabSize = 14, titleLabSize = 20,
                    returnPlot = FALSE) +
  geom_label(aes(20, 80, label = '80% explained variation', 
                 vjust = -1, size = 8))

# Draw multiple bi-plots.
ppairs <- pairsplot(p, components = getComponents(p, c(1:3)),
                    triangle = TRUE, trianglelabSize = 12,
                    hline = 0, vline = 0,
                    pointSize = 0.8, gridlines.major = FALSE, 
                    gridlines.minor = FALSE,
                    colby = 'HealthStatus',
                    title = '', plotaxes = FALSE,
                    margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
                    returnPlot = FALSE)

# Comparing 2 selected principal components
pbiplot <- biplot(p,
                  # loadings parameters
                  showLoadings = TRUE,
                  lengthLoadingsArrowsFactor = 1.5,
                  sizeLoadingsNames = 4,
                  colLoadingsNames = 'red4',
                  # other parameters
                  lab = NULL,
                  colby = 'HealthStatus',  colkey = c('Normal'='royalblue', 
                                                      'Tumor'='red3'),
                  hline = 0, vline = c(-25, 0, 25),
                  vlineType = c('dotdash', 'solid','dashed'),
                  gridlines.major = FALSE, gridlines.minor = FALSE,
                  pointSize = 5,
                  legendPosition = 'none', legendLabSize = 16, 
                  legendIconSize = 8.0,
                  shape = 'Plex',shapekey = 
                    c('01CPTAC_OVprospective_P_PNNL_20161212'=15, 
                    '02CPTAC_OVprospective_P_PNNL_20161212'=17, 
                    '03CPTAC_OVprospective_P_PNNL_20161212'=8, 
                    '04CPTAC_OVprospective_P_PNNL_20161212'=7,
                    '05CPTAC_OVprospective_P_PNNL_20161212'=5, 
                    '06CPTAC_OVprospective_P_PNNL_20161212'=6,
                    '07CPTAC_OVprospective_P_PNNL_20161212'=3, 
                    '08CPTAC_OVprospective_P_PNNL_20161212'=2,
                    '09CPTAC_OVprospective_P_PNNL_20161212'=1, 
                    '10CPTAC_OVprospective_P_PNNL_20161212'=16,
                    '11CPTAC_OVprospective_P_PNNL_20161212'=18, 
                    '12CPTAC_OVprospective_P_PNNL_20161212'=13),
                  drawConnectors = FALSE,
                  title = 'PCA bi-plot',
                  subtitle = 'PC1 versus PC2',
                  caption = '32 PCs ≈ 80%',
                  returnPlot = FALSE)
# Plot the component loadings for selected principal components and label variables driving variation along these.
ploadings <- plotloadings(p, rangeRetain = 0.01, labSize = 4,
                          title = 'Loadings plot', axisLabSize = 12,
                          subtitle = 'PC1, PC2, PC3, PC4, PC5',
                          caption = 'Top 1% variables',
                          shape = 24, shapeSizeRange = c(4, 8),
                          col = c('limegreen', 'black', 'red3'),
                          legendPosition = 'none',
                          drawConnectors = FALSE,
                          returnPlot = FALSE)

# Correlate principal components to continuous variable metadata and test significancies of these.
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c('HealthStatus','Race','Organ',
                                       'TumorStage','Plex'),
                          cexCorval = 1.0,
                          fontCorval = 2,
                          posLab = 'all', 
                          rotLabX = 45,
                          scale = TRUE,
                          main = "PC clinical correlates",
                          cexMain = 1.5,
                          plotRsquared = FALSE,
                          corFUN = 'pearson',
                          corUSE = 'pairwise.complete.obs',
                          signifSymbols = c('****', '***', '**', '*', ''),
                          signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                          returnPlot = FALSE)

library(cowplot)
library(ggplotify)

# Plot everything together
top_row <- plot_grid(pscree, ppairs, pbiplot,
                     ncol = 3,
                     labels = c('A', 'B  Pairs plot', 'C'),
                     label_fontfamily = 'serif',
                     label_fontface = 'bold',
                     label_size = 22,
                     align = 'h',
                     rel_widths = c(1.10, 0.80, 1.10))

bottom_row <- plot_grid(ploadings,
                        as.grob(peigencor),
                        ncol = 2,
                        labels = c('D', 'E'),
                        label_fontfamily = 'serif',
                        label_fontface = 'bold',
                        label_size = 22,
                        align = 'h',
                        rel_widths = c(0.8, 1.2))

plot_grid(top_row, bottom_row, ncol = 1,
          rel_heights = c(1.1, 0.9))


pbiplot <- biplot(p,
                  # loadings parameters
                  showLoadings = TRUE,
                  lengthLoadingsArrowsFactor = 1.5,
                  sizeLoadingsNames = 4,
                  colLoadingsNames = 'red4',
                  # other parameters
                  lab = NULL,
                  colby = 'TumorStage',
                  hline = 0, vline = c(-25, 0, 25),
                  vlineType = c('dotdash', 'solid','dashed'),
                  gridlines.major = FALSE, gridlines.minor = FALSE,
                  pointSize = 5,
                  legendPosition = 'right', legendLabSize = 16, 
                  legendIconSize = 8.0,
                  shape = 'Plex',
                  drawConnectors = FALSE,
                  title = 'PCA bi-plot',
                  subtitle = 'PC1 versus PC2',
                  caption = '27 PCs ≈ 80%',
                  returnPlot = FALSE)


# Check the optimal number of cluster for CPTAC
fviz_nbclust(CPTAC,kmeans, method="silhouette")

# Check for batch effect in the plex information and corrected
batch=design$Plex

# HealthStatus, Race, Organ and TumorStage as predictors in the model
design0=model.matrix(~HealthStatus+Race+Organ+TumorStage, data=design)

# Correct batch effect from CPTAC dataset using batch information and model matrix 
remov_batch=removeBatchEffect(CPTAC,batch,design=design0)

# Check the batch effect with boxplots
par(mfrow=c(1,2))
boxplot(as.data.frame(CPTAC),main="Original")
boxplot(as.data.frame(remov_batch),main="Batch Corrected")

# Boxplots show that ther is no batch effect

# Create distict colors for categories creating annotated heatmaps
generate_annotation_colors <- function(design, columns) {
  annotation_colors <- list()
  for (col_name in columns) {
    n_categories <- length(unique(design[[col_name]]))
    col_palette <- 
      colorRampPalette(grDevices::rainbow(n_categories))(n_categories)
    annotation_colors[[col_name]] <- setNames(col_palette, 
                                              unique(design[[col_name]]))
  }
  return(annotation_colors)
}

# Specify the columns you want to generate annotation colors for
columns_to_color <- c("HealthStatus", "Race", "Organ","TumorStage","Plex")

# Generate annotation colors for the specified columns
mycolors <- generate_annotation_colors(design, columns_to_color)

# Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
hc <- hclust(dist(CPTAC)) # 3 clusters

# Plot hierarchical clustering
plot(hc,hang = -1,labels = FALSE,xlab = "",sub="")

# Create a pheatmap with 3 clusters as sugestinh from hc
pheat_data=pheatmap(CPTAC, annotation_col = design,
                    clustering_distance_cols = "euclidean",
                    cutree_cols = 3,
                    show_rownames = F,annotation_colors=mycolors)

# Name the clusters
cluster=cutree(pheat_data$tree_col,3)
annotation2 = data.frame(cluster)
rownames(annotation2) = colnames(CPTAC)

# Merge annotation files with the name of the clusters
design2=merge(design,annotation2,by="row.names")
design2
rownames(design2) <- design2[,1]
design2=design2[,-1]


# Final heatmap with clusters
pheatmap(CPTAC, annotation_col = design2,clustering_distance_cols = "euclidean",cutree_cols = 3,
         show_rownames = F,annotation_colors=mycolors)

# Make CPTAC dataset as dataframe
CPTAC=as.data.frame(CPTAC)


# Extract samples based on the three clusters
cluster_1_samples <- rownames(design2)[design2$cluster == 1]
cluster_2_samples <- rownames(design2)[design2$cluster == 2]
cluster_3_samples <- rownames(design2)[design2$cluster == 3]

# Preparing a subset of the CPTAC data for cluster 1, reformatting it, and labeling it for identification
gene_df_cluster_1 <- CPTAC[, cluster_1_samples]
gene_df_cluster_1$Gene <- rownames(gene_df_cluster_1)
gene_df_cluster_1 <- gene_df_cluster_1 %>% pivot_longer(-Gene) %>%
  pivot_wider(names_from="Gene", values_from="value") %>%
  rename(id=name)
gene_df_cluster_1<-gene_df_cluster_1%>%
  mutate(group="1")

# Preparing a subset of the CPTAC data for cluster 2, reformatting it, and labeling it for identification
gene_df_cluster_2 <- CPTAC[, cluster_2_samples]
gene_df_cluster_2$Gene <- rownames(gene_df_cluster_2)
gene_df_cluster_2 <- gene_df_cluster_2 %>% pivot_longer(-Gene) %>%
  pivot_wider(names_from="Gene", values_from="value") %>%
  rename(id=name)
gene_df_cluster_2 <-gene_df_cluster_2%>%
  mutate(group="2")

# Preparing a subset of the CPTAC data for cluster 3, reformatting it, and labeling it for identification
gene_df_cluster_3 <- CPTAC[, cluster_3_samples]
gene_df_cluster_3$Gene <- rownames(gene_df_cluster_3)
gene_df_cluster_3 <- gene_df_cluster_3 %>% pivot_longer(-Gene) %>%
  pivot_wider(names_from="Gene", values_from="value") %>%
  rename(id=name)
gene_df_cluster_3<-gene_df_cluster_3%>%
  mutate(group="3")

# Combine datasets by row
merge <- rbind(gene_df_cluster_1, gene_df_cluster_2, gene_df_cluster_3)

#SND1 marker
#One-way ANOVA 
model1=aov(SND1~group,data = merge)
summary(model1)

#Check for homoscedasticity
par(mfrow=c(2,2))
plot(model1)
par(mfrow=c(1,1))

# Post-hoc test to check if there are differences between the group means but not what the differences are.
tukey.model1<-TukeyHSD(model1)
tukey.model1
TK_data_snd1<-as.data.frame(tukey.model1[1])
TK_data_snd1 <- tibble::rownames_to_column(TK_data_snd1, "Group comparison")

# Plot the results in a graph
tukey.plot.aov<-aov(SND1 ~ group, data=merge)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

# Summarize the original data 
mean.SND1.data <- merge %>%
  group_by(group) %>%
  summarise(
    SND1 = mean(SND1)
  )
# Plot the raw data
one.way.plot.snd1 <- ggplot(merge, aes(x = group, y = SND1, )) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))

# Add the means and standard errors to the graph
one.way.plot.snd1 <- one.way.plot.snd1 +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=mean.SND1.data, aes(x=group, y=SND1))
one.way.plot.snd1 <- one.way.plot.snd1 +
  geom_text(data=mean.SND1.data, label=mean.SND1.data$group,
            vjust = -8, size = 5) 
one.way.plot.snd1
one.way.plot.snd1 <- one.way.plot.snd1 +
  theme_classic2() +
  labs(title = "SND1 marker in response to different groups",
       x = "Groups of patients (1=Solid Normal Tissue, 
       2=Tumor Tissue,3=Tumor Tissue)",
       y = "SND2 Gene expression")
one.way.plot.snd1


# MTDH marker
#One-way ANOVA 
model2=aov(MTDH~group,data = merge)
summary(model2)

# Check for homoscedasticity
par(mfrow=c(2,2))
plot(model2)
par(mfrow=c(1,1))

# Post-hoc test to check if there are differences between the group means but not what the differences are.
tukey.model2<-TukeyHSD(model2)
tukey.model2
TK_data_MTDH<-as.data.frame(tukey.model2[1])
TK_data_MTDH <- tibble::rownames_to_column(TK_data_MTDH, "Group comparison")

# Plot the results in a graph
tukey.plot.aov<-aov(MTDH ~ group, data=merge)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

# Summarize the original data 
mean.MTDH.data <- merge %>%
  group_by(group) %>%
  summarise(
    MTDH = mean(MTDH)
  )

#Plot the raw data
one.way.plot.MTDH <- ggplot(merge, aes(x = group, y = MTDH, )) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))

#Add the means and standard errors to the graph
one.way.plot.MTDH <- one.way.plot.MTDH +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=mean.MTDH.data, aes(x=group, y=MTDH))
one.way.plot.MTDH <- one.way.plot.MTDH +
  geom_text(data=mean.MTDH.data, label=mean.MTDH.data$group, 
            vjust = -8, size = 5) 
one.way.plot.MTDH
one.way.plot.MTDH <- one.way.plot.MTDH +
  theme_classic2() +
  labs(title = "MTDH marker in response to different groups",
       x = "Groups of patients (1=Solid Normal Tissue,
       2=Tumor Tissue,3=Tumor Tissue)",
       y = "MTDH Gene expression")
one.way.plot.MTDH

#MKI67 marker
#One-way ANOVA 
model3=aov(MKI67~group,data = merge)
summary(model3)

#Check for homoscedasticity
par(mfrow=c(2,2))
plot(model3)
par(mfrow=c(1,1))

# Post-hoc test to check if there are differences between the group means but not what the differences are.
tukey.model3<-TukeyHSD(model3)
tukey.model3
TK_data_MKI67<-as.data.frame(tukey.model3[1])
TK_data_MKI67 <- tibble::rownames_to_column(TK_data_MKI67, "Group comparison")

# Plot the results in a graph
tukey.plot.aov<-aov(MKI67 ~ group, data=merge)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

# Summarize the original data 
mean.MKI67.data <- merge %>%
  group_by(group) %>%
  summarise(
    MKI67 = mean(MKI67)
  )
# Plot the raw data
one.way.plot.MKI67 <- ggplot(merge, aes(x = group, y = MKI67, )) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))

# Add the means and standard errors to the graph
one.way.plot.MKI67 <- one.way.plot.MKI67 +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=mean.MKI67.data, aes(x=group, y=MKI67))
one.way.plot.MKI67 <- one.way.plot.MKI67 +
  geom_text(data=mean.MKI67.data, label=mean.MKI67.data$group, 
            vjust = -8, size = 5) 
one.way.plot.MKI67
one.way.plot.MKI67 <- one.way.plot.MKI67 +
  theme_classic2() +
  labs(title = "MKI67 marker in response to different groups",
       x = "Groups of patients (1=Solid Normal Tissue, 2=Tumor Tissue,
       3=Tumor Tissue)",
       y = "MKI67 Gene expression")
one.way.plot.MKI67


# SND1 Turkey multiple comparissons
main.title <- "Tukey multiple comparisons of means"
subtitle <- paste0(
  "95% family-wise confidence level",
  " Fit: aov(formula = SND1 ~ group, data = merge).") %>%
  strwrap(width = 80) %>%
  paste(collapse = "\n")
tab=ggtexttable(TK_data_snd1,theme = ttheme("light"),rows=NULL)
tab%>%
  tab_add_title(text = subtitle, face = "plain", size = 10) %>%
  tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) 

# MTDH Turkey multiple comparissons
main.title <- "Tukey multiple comparisons of means"
subtitle <- paste0(
  "95% family-wise confidence level",
  " Fit: aov(formula = MTDH ~ group, data = merge).") %>%
  strwrap(width = 80) %>%
  paste(collapse = "\n")
tab2=ggtexttable(TK_data_MTDH,theme = ttheme("light"),rows=NULL)
tab2%>%
  tab_add_title(text = subtitle, face = "plain", size = 10) %>%
  tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) 

# MKI67 Turkey multiple comparissons
main.title <- "Tukey multiple comparisons of means"
subtitle <- paste0(
  "95% family-wise confidence level",
  " Fit: aov(formula = MKI67 ~ group, data = merge).") %>%
  strwrap(width = 80) %>%
  paste(collapse = "\n")
tab3=ggtexttable(TK_data_MKI67,theme = ttheme("light"),rows=NULL)
tab3%>%
  tab_add_title(text = subtitle, face = "plain", size = 10) %>%
  tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) 

# Plot everything together
ggarrange(one.way.plot.snd1,one.way.plot.MTDH, one.way.plot.MKI67,tab,tab2,tab3,
          ncol = 3, nrow = 3,
          heights = c(1, 0.5, 0.3))

##===============================================================
##                THIS IS THE END OF THE SCRIPT                ==
##===============================================================
