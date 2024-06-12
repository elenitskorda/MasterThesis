#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
#browseVignettes("ComplexHeatmap")
###BiocManager::install("Seurat")###########CHECK ON THIS
#BiocManager::install("DEGreport")

library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(factoextra)
library (limma)
library(ggpubr)

#Import data
peptide <- read.delim("ratio_peptide_MD.tsv", sep = "\t",check.names = FALSE)
plex <- read.delim("experiment_annotation2.txt", sep = "\t",check.names = FALSE)
clinical = read.delim("PDC_clinical_manifest_09042023_135706.tsv", sep = "\t",check.names = FALSE)
annotation <- read.delim("table_design.tsv", sep = "\t",check.names = FALSE)
clinical2 = read.delim("clinical2.txt", sep = "\t",check.names = FALSE)
clinical2 = subset(clinical2, select = c(1, 4))
sample=clinical2[,1]
rownames(clinical2) = make.names(sample, unique=TRUE)
new <- gsub("X", "", rownames(clinical2))
rownames(clinical2)=new
clinical2=clinical2[, -1]






# Delete columns
peptide <- peptide[, !(names(peptide) %in% c("Index", "ProteinID","Peptide","MaxPepProb","ReferenceIntensity"))]

#Remove NAs
peptide=na.omit(peptide)

Gene=peptide[,1]
#Make unique Gene names
rownames(peptide) = make.names(Gene, unique=TRUE)
#Delete column 1
peptide=peptide[,-1]


#Calculate the variance for each gene in a column called rowVar
peptide$row_var = rowVars(as.matrix(peptide[,c(-1)]))

# Check the range of row var
range <- range(peptide$row_var)
row_var=peptide$row_var
#Sort the values in a decreasing order
sorted_values <- sort(row_var, decreasing = TRUE)
#set an index of 50%
split_index <- length(sorted_values) * 0.5
#Choose the top 50 %
top_50_percent <- sorted_values[1:split_index]
keep_rows <- row_var %in% top_50_percent
filtered_df <- peptide[keep_rows, ]

filtered_df=t(filtered_df)



sample=annotation[,1]



clinical = subset(clinical, select = c(2, 8, 12, 14))
colnames(clinical)[1] ="sample"
colnames(clinical) <- c("sample","Race",'TissueOrgan',"TumorStage")



clinical$sample[clinical$sample == '01OV007'] <- '01OV007.1'
clinical$sample[clinical$sample == '01OV017'] <- '01OV017.1'
clinical$sample[clinical$sample == '01OV023'] <- '01OV023.1'

clinical$sample[clinical$sample == '17OV001'] <- '17OV001.1'

clinical$sample[clinical$sample == '01OV039'] <- '01OV039.1'
clinical$sample[clinical$sample == '11OV002'] <- '11OV002.1'

clinical$sample[clinical$sample == '14OV011'] <- '01OV011.1'

new_row <- data.frame(
  sample = c('01OV029.1','15OV001.1','01OV007','01OV017','01OV023','17OV001','01OV039','11OV002','17OV002.1','14OV011'),
  Race = c("No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data"),
  TissueOrgan= c("No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data"),
  TumorStage=c("No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data","No Data"))
clinical <- rbind(clinical, new_row)

predesign=merge(annotation, clinical,by="sample")
design=merge(predesign, plex,by="sample")
design=subset(design,select=c(1:6))
colnames(design) <- c("sample",'HealthStatus',"Race",'Organ',"TumorStage", "Plex")


row.names(design) <- make.names(design$sample, unique=TRUE)
new <- gsub("X", "", rownames(design))
rownames(design)=new
design=design[,-1]

#prepare heatmap

row_names_match <- all(rownames(design) == rownames(filtered_df))
row_names_df1 <- rownames(filtered_df)
row_names_df2 <- rownames(design)
non_matching_row_names <- setdiff(row_names_df1, row_names_df2)
common_row_names <- intersect(row_names_df1, row_names_df2)

# Filter both data frames to keep only the common row names
df1_filtered <- filtered_df[row_names_df1 %in% common_row_names, ]

row_names_df1 <- rownames(df1_filtered)
row_names_df2 <- rownames(design)
non_matching_row_names <- setdiff(row_names_df1, row_names_df2)

if (length(non_matching_row_names) == 0) {
  cat("All row names match between df1 and df2.\n")
} else {
  cat("Row names that do not match between df1 and df2:\n")
  print(non_matching_row_names)
}
#design=data.matrix(design)


as.matrix(design)
df1_filtered <- df1_filtered[rownames(design), ]

df1_filtered=t(df1_filtered)
row.names(design) <- colnames(df1_filtered)

fviz_nbclust(df1_filtered,kmeans, method="silhouette")

batch=design$Plex


#library(cluster)
#gap.stat <- clusGap(df1_filtered, FUNcluster = kmeans, K.max = 15)
#gap.stat
#fviz_gap_stat(gap.stat)
#library(NbClust)
#NbClust(df1_filtered, method = 'complete', index = 'all')$Best.nc


# Elbow method
#fviz_nbclust(df1_filtered, kmeans, method = "wss") +
#  geom_vline(xintercept = 4, linetype = 2)+
#  labs(subtitle = "Elbow method")
# Silhouette method
#fviz_nbclust(df1_filtered, kmeans, method = "silhouette")+
#  labs(subtitle = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
#set.seed(123)
#fviz_nbclust(df1_filtered, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
#  labs(subtitle = "Gap statistic method")

design0=model.matrix(~HealthStatus+Race+Organ+TumorStage, data=design)
remov_batch=removeBatchEffect(df1_filtered,batch,design=design0)
par(mfrow=c(1,2))
boxplot(as.data.frame(df1_filtered),main="Original")
boxplot(as.data.frame(remov_batch),main="Batch Corrected")

#Testing with pca

sample_data_t <- t(df1_filtered)
pca_result <- prcomp(sample_data_t)
pca_plot=fviz_pca_ind(pca_result, geom=c("point"), habillage=design$HealthStatus)
pca_plot

sample_data_t2 <- t(remov_batch)
pca_result2 <- prcomp(sample_data_t2)
pca_plot2=fviz_pca_ind(pca_result2, geom=c("point"), habillage=design$HealthStatus)
pca_plot2
par(mfrow=c(1,2))
pca_plot
pca_plot2
# Use the combined palette for annotation_colors


generate_annotation_colors <- function(design, columns) {
  annotation_colors <- list()
  for (col_name in columns) {
    n_categories <- length(unique(design[[col_name]]))
    col_palette <- colorRampPalette(grDevices::rainbow(n_categories))(n_categories)
    annotation_colors[[col_name]] <- setNames(col_palette, unique(design[[col_name]]))
  }
  return(annotation_colors)
}

# Specify the columns you want to generate annotation colors for
columns_to_color <- c("HealthStatus", "Race", "Organ","TumorStage","Plex")

# Generate annotation colors for the specified columns
mycolors <- generate_annotation_colors(design, columns_to_color)

drows = dist(df1_filtered, method = "minkowski")
dcols = dist(t(df1_filtered), method = "minkowski")

pheatmap(df1_filtered, annotation_col = design,cutree_cols = 3,cutree_rows = 3,
         show_rownames = F,annotation_colors=mycolors)

pheat_data=pheatmap(df1_filtered, annotation_col = design,cutree_cols = 3,
                    show_rownames = F,annotation_colors=mycolors)

cl = cutree(pheat_data$tree_row,3)
ann = data.frame(cl)
rownames(ann) = rownames(df1_filtered)

col=cutree(pheat_data$tree_col,3)
ann2 = data.frame(col)
rownames(ann2) = colnames(df1_filtered)

design2=merge(design,ann2,by="row.names")
design2
rownames(design2) <- design2[,1]
design2=design2[,-1]

pheatmap(df1_filtered, annotation_col = design2,cutree_cols = 3,
         show_rownames = F,annotation_colors=mycolors)

####################################################################
h <- hclust(df1_filtered, method="complete")
plot( as.dendrogram(h) , las=1, main="d=euclidean\nh=complete")

######################################################################
hc <- hclust(dist(df1_filtered))
pheatmap(df1_filtered, annotation_col = design2,clustering_distance_cols = "euclidean",cutree_cols = 3,
         show_rownames = F,annotation_colors=mycolors)
#PCATOOLS
library(PCAtools)


p <- pca(df1_filtered,metadata = design)

eigencorplot(p,components = getComponents(p, 1:10),metavars = c('HealthStatus','Race','Organ','TumorStage','Plex'))

horn <- parallelPCA(df1_filtered)
horn$n
elbow <- findElbowPoint(p$variance)
elbow
library(ggplot2)

screeplot(p,
          components = getComponents(p, 1:30),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

##############################################################################################################

pscree <- screeplot(p, components = getComponents(p, 1:35),
                    hline = 80, vline = 33, axisLabSize = 14, titleLabSize = 20,
                    returnPlot = FALSE) +
  geom_label(aes(20, 80, label = '80% explained variation', vjust = -1, size = 8))

ppairs <- pairsplot(p, components = getComponents(p, c(1:3)),
                    triangle = TRUE, trianglelabSize = 12,
                    hline = 0, vline = 0,
                    pointSize = 0.8, gridlines.major = FALSE, gridlines.minor = FALSE,
                    colby = 'HealthStatus',
                    title = '', plotaxes = FALSE,
                    margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
                    returnPlot = FALSE)

pbiplot <- biplot(p,
                  # loadings parameters
                  showLoadings = TRUE,
                  lengthLoadingsArrowsFactor = 1.5,
                  sizeLoadingsNames = 4,
                  colLoadingsNames = 'red4',
                  # other parameters
                  lab = NULL,
                  colby = 'HealthStatus',  colkey = c('Normal'='royalblue', 'Tumor'='red3'),
                  hline = 0, vline = c(-25, 0, 25),
                  vlineType = c('dotdash', 'solid','dashed'),
                  gridlines.major = FALSE, gridlines.minor = FALSE,
                  pointSize = 5,
                  legendPosition = 'none', legendLabSize = 16, legendIconSize = 8.0,
                  shape = 'Plex',shapekey = c('01CPTAC_OVprospective_P_PNNL_20161212'=15, 
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

ploadings <- plotloadings(p, rangeRetain = 0.01, labSize = 4,
                          title = 'Loadings plot', axisLabSize = 12,
                          subtitle = 'PC1, PC2, PC3, PC4, PC5',
                          caption = 'Top 1% variables',
                          shape = 24, shapeSizeRange = c(4, 8),
                          col = c('limegreen', 'black', 'red3'),
                          legendPosition = 'none',
                          drawConnectors = FALSE,
                          returnPlot = FALSE)

peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c('HealthStatus','Race','Organ','TumorStage','Plex'),
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
                  legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
                  shape = 'Plex',
                  drawConnectors = FALSE,
                  title = 'PCA bi-plot',
                  subtitle = 'PC1 versus PC2',
                  caption = '27 PCs ≈ 80%',
                  returnPlot = FALSE)


df1_filtered=as.data.frame(df1_filtered)
#install.packages('ggplot2')
#install.packages('GGally')

# #load libraries
# library(ggplot2)
# library(GGally)
# 
# #create pairs plot
# ggpairs(df1_filtered)

#Extract samples based on the clusters.###########################################################
cluster_1_samples <- rownames(design2)[design2$col == 1]
cluster_2_samples <- rownames(design2)[design2$col == 2]
cluster_3_samples <- rownames(design2)[design2$col == 3]

gene_df_cluster_1 <- df1_filtered[, cluster_1_samples]
gene_df_cluster_1$Gene <- rownames(gene_df_cluster_1)
gene_df_cluster_1 <- gene_df_cluster_1 %>% pivot_longer(-Gene) %>%
  pivot_wider(names_from="Gene", values_from="value") %>%
  rename(id=name)
gene_df_cluster_1<-gene_df_cluster_1%>%
  mutate(group="1")

gene_df_cluster_2 <- df1_filtered[, cluster_2_samples]
gene_df_cluster_2$Gene <- rownames(gene_df_cluster_2)
gene_df_cluster_2 <- gene_df_cluster_2 %>% pivot_longer(-Gene) %>%
  pivot_wider(names_from="Gene", values_from="value") %>%
  rename(id=name)
gene_df_cluster_2 <-gene_df_cluster_2%>%
  mutate(group="2")

gene_df_cluster_3 <- df1_filtered[, cluster_3_samples]
gene_df_cluster_3$Gene <- rownames(gene_df_cluster_3)
gene_df_cluster_3 <- gene_df_cluster_3 %>% pivot_longer(-Gene) %>%
  pivot_wider(names_from="Gene", values_from="value") %>%
  rename(id=name)
gene_df_cluster_3<-gene_df_cluster_3%>%
  mutate(group="3")

#merge datasets
merge <- rbind(gene_df_cluster_1, gene_df_cluster_2, gene_df_cluster_3)

#SND1 marker


#One-way ANOVA 
model1=aov(SND1~group,data = merge)
summary(model1)

#Check for homoscedasticity
par(mfrow=c(2,2))
plot(model1)
par(mfrow=c(1,1))
# post-hoc test to check if there are differences between the group means but not what the differences are.

tukey.model1<-TukeyHSD(model1)
tukey.model1
TK_data_snd1<-as.data.frame(tukey.model1[1])
TK_data_snd1 <- tibble::rownames_to_column(TK_data_snd1, "Group comparison")

# Plot the results in a graph
tukey.plot.aov<-aov(SND1 ~ group, data=merge)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

#summarize the original data 
mean.SND1.data <- merge %>%
  group_by(group) %>%
  summarise(
    SND1 = mean(SND1)
  )
#Plot the raw data
one.way.plot.snd1 <- ggplot(merge, aes(x = group, y = SND1, )) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))

#Add the means and standard errors to the graph
one.way.plot.snd1 <- one.way.plot.snd1 +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=mean.SND1.data, aes(x=group, y=SND1))
one.way.plot.snd1 <- one.way.plot.snd1 +
  geom_text(data=mean.SND1.data, label=mean.SND1.data$group, vjust = -8, size = 5) 
one.way.plot.snd1
one.way.plot.snd1 <- one.way.plot.snd1 +
  theme_classic2() +
  labs(title = "SND1 marker in response to different groups",
       x = "Groups of patients (1=Solid Normal Tissue, 2=Tumor Tissue,3=Tumor Tissue)",
       y = "SND2 Gene expression")
one.way.plot.snd1


#MTDH marker


#One-way ANOVA 
model2=aov(MTDH~group,data = merge)
summary(model2)

#Check for homoscedasticity
par(mfrow=c(2,2))
plot(model2)
par(mfrow=c(1,1))
# post-hoc test to check if there are differences between the group means but not what the differences are.

tukey.model2<-TukeyHSD(model2)
tukey.model2
TK_data_MTDH<-as.data.frame(tukey.model2[1])
TK_data_MTDH <- tibble::rownames_to_column(TK_data_MTDH, "Group comparison")
# Plot the results in a graph
tukey.plot.aov<-aov(MTDH ~ group, data=merge)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

#summarize the original data 
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
  geom_text(data=mean.MTDH.data, label=mean.MTDH.data$group, vjust = -8, size = 5) 
one.way.plot.MTDH
one.way.plot.MTDH <- one.way.plot.MTDH +
  theme_classic2() +
  labs(title = "MTDH marker in response to different groups",
       x = "Groups of patients (1=Solid Normal Tissue, 2=Tumor Tissue,3=Tumor Tissue)",
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
# post-hoc test to check if there are differences between the group means but not what the differences are.

tukey.model3<-TukeyHSD(model3)
tukey.model3
TK_data_MKI67<-as.data.frame(tukey.model3[1])
TK_data_MKI67 <- tibble::rownames_to_column(TK_data_MKI67, "Group comparison")
# Plot the results in a graph
tukey.plot.aov<-aov(MKI67 ~ group, data=merge)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

#summarize the original data 
mean.MKI67.data <- merge %>%
  group_by(group) %>%
  summarise(
    MKI67 = mean(MKI67)
  )
#Plot the raw data
one.way.plot.MKI67 <- ggplot(merge, aes(x = group, y = MKI67, )) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))

#Add the means and standard errors to the graph
one.way.plot.MKI67 <- one.way.plot.MKI67 +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=mean.MKI67.data, aes(x=group, y=MKI67))
one.way.plot.MKI67 <- one.way.plot.MKI67 +
  geom_text(data=mean.MKI67.data, label=mean.MKI67.data$group, vjust = -8, size = 5) 
one.way.plot.MKI67
one.way.plot.MKI67 <- one.way.plot.MKI67 +
  theme_classic2() +
  labs(title = "MKI67 marker in response to different groups",
       x = "Groups of patients (1=Solid Normal Tissue, 2=Tumor Tissue,3=Tumor Tissue)",
       y = "MKI67 Gene expression")
one.way.plot.MKI67





#TK_plot_snd1=ggtexttable(TK_data_snd1,theme = ttheme("light"),rows=NULL)
#snd1 Turkey multiple comparissons
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

#MTDH Turkey multiple comparissons
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

#MKI67 Turkey multiple comparissons
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


ggarrange(one.way.plot.snd1,one.way.plot.MTDH, one.way.plot.MKI67,tab,tab2,tab3,
          ncol = 3, nrow = 3,
          heights = c(1, 0.5, 0.3))
#cluste1
design3=design2

sample_names_gene_df <- colnames(gene_df_cluster_1)
design3_filtered <- design3[row.names(design3) %in% sample_names_gene_df, , drop = FALSE]


#cluster2
design4=design2

sample_names_gene_df2 <- colnames(gene_df_cluster_2)
design4_filtered <- design3[row.names(design4) %in% sample_names_gene_df2, , drop = FALSE]

#cluster3

design5=design2

sample_names_gene_df3 <- colnames(gene_df_cluster_3)
design5_filtered <- design3[row.names(design5) %in% sample_names_gene_df3, , drop = FALSE]

###
#differential expression
###

#cluster1
group <- design3_filtered$HealthStatus
# Create a design matrix
design <- model.matrix(~0 + group)

cont_matrix <- makeContrasts(TumorvsNormal = groupTumor-groupNormal,
                             levels=design)
print(cont_matrix)
# Fit the expression matrix to a linear model
fit <- lmFit(gene_df_cluster_1, design)
# Compute contrast
fit_contrast <- contrasts.fit(fit, cont_matrix)
# Bayes statistics of differential expression
# *There are several options to tweak!*
fit_contrast <- eBayes(fit_contrast)
# Get the names of the genes
gene_names <- rownames(fit_contrast$coefficients)
# Generate the volcano plot with gene names as labels
#it highlights the gignificant ones!

#cluster2
group2 <- design4_filtered$HealthStatus
# Create a design matrix
design2 <- model.matrix(~0 + group2)

cont_matrix2 <- makeContrasts(TumorvsNormal = group2Tumor-group2Normal,
                              levels=design2)
print(cont_matrix2)
# Fit the expression matrix to a linear model
fit2 <- lmFit(gene_df_cluster_2, design2)
# Compute contrast
fit_contrast2 <- contrasts.fit(fit2, cont_matrix2)
# Bayes statistics of differential expression
# *There are several options to tweak!*
fit_contrast2 <- eBayes(fit_contrast2)
# Get the names of the genes
gene_names2 <- rownames(fit_contrast2$coefficients)
# Generate the volcano plot with gene names as labels
#it highlights the gignificant ones



##########################################################################


#cluster3
group3 <- design5_filtered$HealthStatus
# Create a design matrix
design3 <- model.matrix(~0 + group3)

cont_matrix3 <- makeContrasts(TumorvsNormal = group3Tumor-group3Normal,
                              levels=design3)
print(cont_matrix2)
# Fit the expression matrix to a linear model
fit2 <- lmFit(gene_df_cluster_2, design2)
# Compute contrast
fit_contrast2 <- contrasts.fit(fit2, cont_matrix2)
# Bayes statistics of differential expression
# *There are several options to tweak!*
fit_contrast2 <- eBayes(fit_contrast2)
# Get the names of the genes
gene_names2 <- rownames(fit_contrast2$coefficients)
# Generate the volcano plot with gene names as labels
#it highlights the gignificant ones

par(mfcol = c(1,2))
volcanoplot(fit_contrast, highlight = 100, names = gene_names)+
  abline(v = 0, col = "gray")

volcanoplot(fit_contrast2, highlight = 100, names = gene_names)+
  abline(v = 0, col = "gray")

volcanoplot(fit_contrast3, highlight = 100, names = gene_names)+
  abline(v = 0, col = "gray")
