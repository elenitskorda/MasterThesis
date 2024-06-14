# -*- coding: utf-8 -*-

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Date: Mon Sep 5 10:28:18 2023                                                                                  ::
##  File Name: proteome_integration.R                                                                              ::
##  Author: Eleni Theofania Skorda                                                                                 ::
##                                                                                                                 ::
##  Description:                                                                                                   ::
##      Integrates proteomic local cell line data (Reference) with CPTAC data (Query)                              ::
##      in order to indetify cell lines within the query dataset.                                                  ::
##                                                                                                                 ::
##  List of Functions:                                                                                             ::
##      none                                                                                                       ::
##                                                                                                                 ::
##  Procedure:                                                                                                     ::
##      1. Prepare reference and query datasets to a proper format                                                 ::
##      2. Performing integration                                                                                  ::
##      3. Finds markers (differentially expressed genes) for each of the identity classes in a dataset            ::
##      4. Performs enrichment analysis for each cluster and find common pathways between benign and no benign     ::
##      5. Find common genes betweenthe shared pathways (benigh and no benign) and all markers from each cluster   ::
##                                                                                                                 ::
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')
# install.packages('BiocManager')
# BiocManager::install('multtest')
# install.packages('metap')


# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(writexl)
library(tidyverse)
library(patchwork)
library(enrichR)
library(presto)
library(multtest)
library(metap)
library(VennDetail)
library(VennDiagram)
library(ggtext)

message <- "Preparation of datasets."

system2(command = "PowerShell", 
        args = c("-Command", 
                 "\"Add-Type -AssemblyName System.Speech;",
                 "$speak = New-Object System.Speech.Synthesis.SpeechSynthesizer;",
                 paste0("$speak.Speak('", message, "');\"")
        ))


##==========================================================================
##  SECTION 1: Prepare reference and query datasets to a proper format    ==
##==========================================================================


#                            REFERENCE DATASET
##---------------------------------------------------------------------------

# Load the local cell line local.cell.lines into R
local.cell.lines <- read.delim("OvCa Proteins 5sub Loess DE_stats.tsv", 
                               sep = "\t", check.names = FALSE)

# Remove unnecessary columns from the local.cell.lines
local.cell.lines = subset(local.cell.lines, select = -c(1:3))
local.cell.lines = subset(local.cell.lines, select = -c(2:12))

# Make Genes as rownames
Gene = local.cell.lines[, 1]
rownames(local.cell.lines) = make.names(Gene, unique = TRUE)
local.cell.lines = local.cell.lines[, -1]
local.cell.lines = as.data.frame(local.cell.lines)

# Transform data back to Raw values
local.cell.lines=2^local.cell.lines

## Impute low value for missing data in local cell lines
# Identify numeric columns
numeric_cols <- sapply(local.cell.lines, is.numeric)

# Apply low value imputation to numeric columns
for(i in which(numeric_cols)) {
  #  # Calculate the minimum value in the column
  min_val <- min(local.cell.lines[, i], na.rm = TRUE)
  
  # Define a low value
  low_val <- min_val
  
  #  # Replace NA values with the low value
  local.cell.lines[is.na(local.cell.lines[, i]), i] <- low_val
}
##---------------------------------------------------------------------------
#                              QUERY DATASET
##---------------------------------------------------------------------------

# Read query data and remove unnecessary columns
CPTAC <- read.delim("abundance_gene_MD_18.tsv", 
                      sep = "\t", check.names = FALSE)
CPTAC <- CPTAC[, !(names(CPTAC) %in% c("NumberPSM", 
                                             "ProteinID",
                                             "CPTAC", 
                                             "MaxPepProb", 
                                             "ReferenceIntensity"))]

# Make unique Gene names
Gene = CPTAC[, 1]
rownames(CPTAC) = make.names(Gene, unique = TRUE)
CPTAC = CPTAC[, -1]

# Transform data back to Raw values
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

# Read annotation file for reference data
ref.annotation <- read.delim("annotation.txt", sep = "\t", check.names = FALSE)
ref.annotation = ref.annotation %>% remove_rownames %>% 
  column_to_rownames(var = "Sample")

# Read annotation file for quey data
query.annotation <- read.delim("table_design.tsv", sep = "\t", 
                               check.names = FALSE)
query.annotation$group <- gsub("_", "", query.annotation$group)
query.annotation$new_sample <- paste(query.annotation$group, 
                                     query.annotation$sample, sep = "_")

# Create a mapping between old and new column names
column_mapping <- setNames(query.annotation$new_sample, 
                           query.annotation$sample)

# Rename the columns in the query dataset
colnames(CPTAC) <- column_mapping[colnames(CPTAC)]

query.annotation = query.annotation[, -1]
query.annotation = query.annotation %>% 
  remove_rownames %>% 
  column_to_rownames(var = "new_sample")

##---------------------------------------------------------------------------
message <- "Integration."

system2(command = "PowerShell", 
        args = c("-Command", 
                 "\"Add-Type -AssemblyName System.Speech;",
                 "$speak = New-Object System.Speech.Synthesis.SpeechSynthesizer;",
                 paste0("$speak.Speak('", message, "');\"")
        ))

##================================================================
##                   SECTION 2: Integration                     ==
##================================================================


# Create Seurat object for reference data
reference = CreateSeuratObject(
  counts = local.cell.lines,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = ref.annotation,
  project = "my_project"
)

# Create Seurat object for query data
query = CreateSeuratObject(
  counts = CPTAC,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = query.annotation,
  project = "my_project2"
)

# Pre-process the reference dataset without integration
reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference, approx = FALSE)
reference <- FindNeighbors(reference, dims = 1:4)
reference <- FindClusters(reference)
head(Idents(reference), 54)
reference <- RunUMAP(reference, dims = 1:5,reduction = "pca")
DimPlot(reference, group.by = c("Cell_type"), label = TRUE)

# Normalize query dataset
query <- NormalizeData(query)

# Find anchors for integration
anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:6, 
                               reference.reduction = "pca",k.anchor=10)
predictions <- TransferData(anchorset = anchors, refdata = reference$Cell_type,
                            dims = 1:6, k.weight = 6)

# Add predicted cell types to the query dataset
ovarian.query <- AddMetaData(query, metadata = predictions)

# Evaluate the quality of the predictions
ovarian.query$prediction.match <- ovarian.query$predicted.id == ovarian.query$group
table(ovarian.query$prediction.match)

# Display the distribution of predicted cell types
table(ovarian.query$predicted.id)

# Visualize expression of specific genes in predicted cell types
VlnPlot(ovarian.query, c("MTDH", "SND1"), group.by = "predicted.id", pt.size=1.5)

# Run UMAP on the reference dataset for visualization
ovarian.ref <- RunUMAP(reference, dims = 1:6, return.model = TRUE)

# Map the query dataset to the reference UMAP
ovarian.query <- MapQuery(anchorset = anchors, 
                          reference = ovarian.ref, query = ovarian.query, 
                          refdata = list(predicted_celltype = "Cell_type"),
                          reference.reduction = "pca", reduction.model = "umap", 
                          transferdata.args = list(k.weight = 6))

# Plot side-by-side UMAPs of the reference and query datasets
p1 <- DimPlot(ovarian.ref, reduction = "umap", group.by = "Cell_type", 
              label = TRUE, label.size = 5, pt.size = 5, repel = TRUE) + 
              ggtitle("Reference annotations")
p2 <- DimPlot(ovarian.query, reduction = "ref.umap",
              group.by = "predicted.predicted_celltype", label = TRUE,
              label.size = 5, pt.size = 3, repel = TRUE) + 
              ggtitle("Query transferred labels")
p3 <- DimPlot(ovarian.query,pt.size =1)+ggtitle("Query annotations")
plot <- p1+p3+p2+plot_annotation(title = "Proteome_VS_Proteome")  &  
              theme(plot.title = element_text(hjust = 0.5,size = 15)) 
plot

message <- "Exploration of markers"

system2(command = "PowerShell", 
        args = c("-Command", 
                 "\"Add-Type -AssemblyName System.Speech;",
                 "$speak = New-Object System.Speech.Synthesis.SpeechSynthesizer;",
                 paste0("$speak.Speak('", message, "');\"")
        ))

##=============================================================================================================
##  SECTION 3: Finds markers (differentially expressed genes) for each of the identity classes in a dataset  ==
##=============================================================================================================

# Rename the idents of ovarian.quert with integrated ones
Idents(ovarian.query)<-ovarian.query$predicted.id
Idents(ovarian.query)

# Exclude SKOV3 cell line and create ovarian.tumors as a no benign Seurat Obj
ovarian.tumors= subset(x = (ovarian.query), idents = 'SKOV3', invert = TRUE)

# Merge Reference and query
# Keep only certain aspects of the Seurat object. Useful for merge()
reference <- DietSeurat(ovarian.ref, counts = FALSE, data = TRUE,
                       scale.data = FALSE, dimreducs = "umap")
Idents(reference)=reference$Cell_type

query <- DietSeurat(ovarian.query, counts = FALSE, data = TRUE,
                   scale.data = FALSE, dimreducs = "ref.umap", assays="RNA")
reference$id <- 'reference'
query$id <- 'query'

# Merge the reference and query Seurat Obj
reference.query <- merge(reference, y=query)

# Merge the reductions
reference.query[["umap"]] <- merge(reference[["umap"]], query[["ref.umap"]])

DefaultAssay(reference.query) <- "RNA"
# Name the ID of the merged dataset as "groups"
reference.query[['groups']] <-reference.query$id

# Split and Join Layers Together
reference.query=JoinLayers(reference.query)

# Finds differentially expressed genes for each of the identity classes in a dataset
all.markers <- FindAllMarkers(object = ovarian.query,only.pos = FALSE)

# Exported it to a csv file for later use
write.csv(all.markers, file="all_markers_proteome.csv")

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
P1

message <- "Enrichment analysis."

system2(command = "PowerShell", 
        args = c("-Command", 
                 "\"Add-Type -AssemblyName System.Speech;",
                 "$speak = New-Object System.Speech.Synthesis.SpeechSynthesizer;",
                 paste0("$speak.Speak('", message, "');\"")
        ))

##===================================================================
##    SECTION 4: Performs enrichment analysis for each cluster     ==
##         and find common pathways between benign and no benign   ==
##===================================================================


##----------------------------------------------------------------
##                   KURAM Cell Line - BENIGN                   --
##----------------------------------------------------------------

# DE and EnrichR pathway return.gene.list = TRUE, no visualization for reference
ref=DEenrichRPlot(reference,
                  ident.1 = "KURAM",
                  max.genes = 100,
                  test.use = "wilcox",
                  p.val.cutoff = 0.05,
                  cols = NULL,
                  enrich.database = "MSigDB_Oncogenic_Signatures",
                  num.pathway = 30,
                  logfc.threshold = 0.1,
                  return.gene.list = T)

# DE and EnrichR pathway return.gene.list = TRUE, no visualization for query
query=DEenrichRPlot(ovarian.query, 
                    ident.1 = "KURAM",
                    max.genes = 100,
                    test.use = "wilcox",
                    p.val.cutoff = 0.05,
                    cols = NULL,
                    enrich.database = "MSigDB_Oncogenic_Signatures",
                    num.pathway = 30,
                    logfc.threshold = 0.1,
                    return.gene.list = T)


# Extract positive markers for reference and query
venpos<- venndetail(list(reference =ref$pos$MSigDB_Oncogenic_Signatures.Term ,
                         Query = query$pos$MSigDB_Oncogenic_Signatures.Term
))

# Extract negative markers for reference and query
venneg<- venndetail(list(reference =ref$neg$MSigDB_Oncogenic_Signatures.Term ,
                         Query = query$neg$MSigDB_Oncogenic_Signatures.Term
))

#Find shared pathways between reference and query for benign data
allpos=getSet(venpos, subset = c("Shared"))

# Create columns and fill them up to create a structured dataframe
allpos %>% 
  mutate(ID="Pos") %>% 
  mutate(Status="BENIGN") %>% 
  mutate(CellLine="KURAM") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->allposKURAM

#Find shared pathways between reference and query for benign data
allneg=getSet(venneg, subset = c("Shared"))  

# Create columns and fill them up to create a structured dataframe
allneg %>% 
  mutate(ID="Neg") %>% 
  mutate(Status="BENIGN") %>% 
  mutate(CellLine="KURAM") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->allnegKURAM

# Combine datasets by row
KURAM_BENIGN <- rbind(allposKURAM, allnegKURAM)

# Extract genes for comparison to expression analysis
query$pos %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, 
           MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryposgenes

allposKURAM %>% 
  left_join(queryposgenes, by="Pathway")-> KURAM_pos_benign

KURAM_pos_benign

query$neg %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, 
           MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryneggenes

allnegKURAM %>% 
  left_join(queryneggenes, by="Pathway")-> KURAM_neg_benign

KURAM_GENES_BENIGN<- rbind(KURAM_neg_benign, KURAM_pos_benign)

##---------------------------------------------------------------
##                 KURAM Cell Line - NO BENIGN                 --
##---------------------------------------------------------------


ref=DEenrichRPlot(reference,
                  ident.1 = "KURAM",
                  max.genes = 100,
                  test.use = "wilcox",
                  p.val.cutoff = 0.05,
                  cols = NULL,
                  enrich.database = "MSigDB_Oncogenic_Signatures",
                  num.pathway = 30,
                  logfc.threshold = 0.1,
                  return.gene.list = T)

query.tumors=DEenrichRPlot(ovarian.tumors, 
                           ident.1 = "KURAM",
                           max.genes = 100,
                           test.use = "wilcox",
                           p.val.cutoff = 0.05,
                           cols = NULL,
                           enrich.database = "MSigDB_Oncogenic_Signatures",
                           num.pathway = 30,
                           logfc.threshold = 0.1,
                           return.gene.list = T)

ven.tumors.pos <- venndetail(list(reference = 
                                  ref$pos$MSigDB_Oncogenic_Signatures.Term , 
                                  Query = 
                                  query.tumors$pos$MSigDB_Oncogenic_Signatures.Term
))
ven.tumors.neg <- venndetail(list(reference =
                                  ref$neg$MSigDB_Oncogenic_Signatures.Term ,
                                  Query = 
                                  query.tumors$neg$MSigDB_Oncogenic_Signatures.Term
))
#Do the same for Tumour data
tumors.pos=getSet(ven.tumors.pos, subset = c("Shared"))

tumors.pos %>% 
  mutate(ID="Pos") %>% 
  mutate(Status="NO_BENIGN") %>% 
  mutate(CellLine="KURAM") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->tumors.posKURAM

tumors.neg=getSet(ven.tumors.neg, subset = c("Shared"))

tumors.neg %>% 
  mutate(ID="Neg") %>% 
  mutate(Status="NO_BENIGN") %>% 
  mutate(CellLine="KURAM") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->tumors.negKURAM

KURAM_NOBENIGN <- rbind(tumors.posKURAM, tumors.negKURAM)
KURAM=rbind(KURAM_BENIGN,KURAM_NOBENIGN)

#Extract genes as well for comparison to expression analysis
query.tumors$pos %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryposgenes

tumors.posKURAM %>% 
  left_join(queryposgenes, by="Pathway") -> KURAM_pos_nobenign

query.tumors$neg %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryneggenes

tumors.negKURAM %>% 
  left_join(queryneggenes, by="Pathway")-> KURAM_neg_nobenign

KURAM_GENES_NOBENIGN<- rbind(KURAM_neg_nobenign, KURAM_pos_nobenign)

#Make one dataframe with genes for KURAM
KURAM_GENES<-rbind(KURAM_GENES_BENIGN, KURAM_GENES_NOBENIGN)

#Remove duplicate pathway entries for both neg and pos
KURAM %>%
  group_by(Pathway, CellLine, Status) %>% 
  count() %>% 
  filter(n>1)
KURAM_alt <- KURAM[!(KURAM$Pathway == "RB DN.V1 DN" & KURAM$ID == "Neg"), ]

KURAM_alt

KURAM_alt %>% 
  pivot_wider(names_from=Status, values_from=ID) %>% 
  mutate(Pathway_bold=ifelse(BENIGN==NO_BENIGN, paste("**", Pathway, "**", sep=""), Pathway)) %>% 
  mutate(Pathway_bold=ifelse(is.na(Pathway_bold), Pathway, Pathway_bold)) %>% 
  mutate(Pathway=Pathway_bold) %>% 
  select(!Pathway_bold) %>% 
  pivot_longer(c(BENIGN, NO_BENIGN), names_to="Status", values_to="ID") %>% 
  drop_na()-> KURAM_plot

KURAM_plot

p1 <- ggplot(KURAM_plot, aes(x=Status , y=Pathway, colour = ID)) +
  geom_point(alpha=2,size = 6)+
  theme_classic()+
  theme(axis.text.y = element_markdown(),plot.title = element_text(hjust = 0.5),
        legend.position = "none", axis.title.x=element_blank())+
        ggtitle("KURAM")

p1


##----------------------------------------------------------------
##                   CAOV3 Cell Line - BENIGN                   --
##----------------------------------------------------------------

# DE and EnrichR pathway return.gene.list = TRUE, no visualization for reference
ref=DEenrichRPlot(reference,
                  ident.1 = "CAOV3",
                  max.genes = 100,
                  test.use = "wilcox",
                  p.val.cutoff = 0.05,
                  cols = NULL,
                  enrich.database = "MSigDB_Oncogenic_Signatures",
                  num.pathway = 30,
                  logfc.threshold = 0.1,
                  return.gene.list = T)
# DE and EnrichR pathway return.gene.list = TRUE, no visualization for query
query=DEenrichRPlot(ovarian.query, 
                    ident.1 = "CAOV3",
                    max.genes = 100,
                    test.use = "wilcox",
                    p.val.cutoff = 0.05,
                    cols = NULL,
                    enrich.database = "MSigDB_Oncogenic_Signatures",
                    num.pathway = 30,
                    logfc.threshold = 0.1,
                    return.gene.list = T)

# Extract positive markers for reference and query
venpos<- venndetail(list(reference =ref$pos$MSigDB_Oncogenic_Signatures.Term , 
                         Query = query$pos$MSigDB_Oncogenic_Signatures.Term
))

# Extract negative markers for reference and query
venneg<- venndetail(list(reference =ref$neg$MSigDB_Oncogenic_Signatures.Term , 
                         Query = query$neg$MSigDB_Oncogenic_Signatures.Term
))

# Extract the common positive marker between reference and query
allpos=getSet(venpos, subset = c("Shared"))

# Create columns and fill them up to create a structured dataframe
allpos %>% 
  mutate(ID="Pos") %>% 
  mutate(Status="BENIGN") %>% 
  mutate(CellLine="CAOV3") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->allposCAOV3

# Extract the common negative marker between reference and query
allneg=getSet(venneg, subset = c("Shared"))  

# Create columns and fill them up to create a structured dataframe
allneg %>% 
  mutate(ID="Neg") %>% 
  mutate(Status="BENIGN") %>% 
  mutate(CellLine="CAOV3") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->allnegCAOV3

# Combine datasets by row
CAOV3_BENIGN <- rbind(allposCAOV3, allnegCAOV3)

# Extract genes as well for comparison to expression analysis
# Both positive and negative
# Positive
query$pos %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, 
           MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryposgenes

allposCAOV3 %>% 
  left_join(queryposgenes, by="Pathway")-> CAOV3_pos_benign

# Negative
query$neg %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, 
           MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryneggenes

allnegCAOV3 %>% 
  left_join(queryneggenes, by="Pathway")-> CAOV3_neg_benign

# Combine datasets by row can create CAOV3_GENES_BENIGN
CAOV3_GENES_BENIGN<- rbind(CAOV3_neg_benign, CAOV3_pos_benign)

##---------------------------------------------------------------
##                 CAOV3 Cell Line - NO BENIGN                 --
##---------------------------------------------------------------

ref=DEenrichRPlot(reference,
                  ident.1 = "CAOV3",
                  max.genes = 100,
                  test.use = "wilcox",
                  p.val.cutoff = 0.05,
                  cols = NULL,
                  enrich.database = "MSigDB_Oncogenic_Signatures",
                  num.pathway = 30,
                  logfc.threshold = 0.1,
                  return.gene.list = T)

query.tumors=DEenrichRPlot(ovarian.tumors, 
                           ident.1 = "CAOV3",
                           max.genes = 100,
                           test.use = "wilcox",
                           p.val.cutoff = 0.05,
                           cols = NULL,
                           enrich.database = "MSigDB_Oncogenic_Signatures",
                           num.pathway = 30,
                           logfc.threshold = 0.1,
                           return.gene.list = T)

ven.tumors.pos <- venndetail(list(reference =
                                  ref$pos$MSigDB_Oncogenic_Signatures.Term , 
                                  Query = 
                                  query.tumors$pos$MSigDB_Oncogenic_Signatures.Term
))
ven.tumors.neg <- venndetail(list(reference =
                                  ref$neg$MSigDB_Oncogenic_Signatures.Term ,
                                  Query = 
                                  query.tumors$neg$MSigDB_Oncogenic_Signatures.Term
))

#Do the same for Tumour data
tumors.pos=getSet(ven.tumors.pos, subset = c("Shared"))

tumors.pos %>% 
  mutate(ID="Pos") %>% 
  mutate(Status="NO_BENIGN") %>% 
  mutate(CellLine="CAOV3") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->tumors.posCAOV3

tumors.neg=getSet(ven.tumors.neg, subset = c("Shared"))

tumors.neg %>% 
  mutate(ID="Neg") %>% 
  mutate(Status="NO_BENIGN") %>% 
  mutate(CellLine="CAOV3") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->tumors.negCAOV3

CAOV3_NOBENIGN <- rbind(tumors.posCAOV3, tumors.negCAOV3)
CAOV3=rbind(CAOV3_BENIGN, CAOV3_NOBENIGN)

#Extract genes as well for comparison to expression analysis
query.tumors$pos %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryposgenes

tumors.posCAOV3 %>% 
  left_join(queryposgenes, by="Pathway")-> CAOV3_pos_nobenign

query.tumors$neg %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryneggenes

tumors.negCAOV3 %>% 
  left_join(queryneggenes, by="Pathway")-> CAOV3_neg_nobenign

CAOV3_GENES_NOBENIGN<- rbind(CAOV3_neg_nobenign, CAOV3_pos_nobenign)
CAOV3_GENES<-rbind(CAOV3_GENES_NOBENIGN, CAOV3_GENES_BENIGN)

#Check which ones have conflicting results
CAOV3 %>%
  group_by(Pathway, CellLine, Status) %>% 
  count() %>% 
  filter(n>1)
CAOV3 <- CAOV3[!(CAOV3$Pathway == "MEK UP.V1 UP" & CAOV3$ID == "Pos"), ]
CAOV3 <- CAOV3[!(CAOV3$Pathway == "RB DN.V1 DN" & 
                   CAOV3$Status=="BENIGN" & CAOV3$ID == "Neg"), ]
CAOV3_alt <- CAOV3[!(CAOV3$Pathway == "RAF UP.V1 UP" & CAOV3$ID == "Pos"), ]

CAOV3_alt
CAOV3_alt %>% 
  pivot_wider(names_from=Status, values_from=ID) %>% 
  mutate(Pathway_bold=ifelse(BENIGN==NO_BENIGN, 
                             paste("**", Pathway, "**", sep=""), Pathway)) %>% 
  mutate(Pathway_bold=ifelse(is.na(Pathway_bold), Pathway, Pathway_bold)) %>% 
  mutate(Pathway=Pathway_bold) %>% 
  select(!Pathway_bold) %>% 
  pivot_longer(c(BENIGN, NO_BENIGN), names_to="Status", values_to="ID") %>% 
  drop_na()-> CAOV3_plot

CAOV3_plot

p2 <- ggplot(CAOV3_plot, aes(x=Status , y=Pathway, colour = ID)) +
  geom_point(alpha=2,size = 6)+
  theme_classic()+
  theme(axis.text.y = element_markdown(),plot.title = element_text(hjust = 0.5),
        legend.position = "none", axis.title.y=element_blank())+
  ggtitle("CAOV3")

p2

##---------------------------------------------------------------
##                  OVSAHO Cell Line - BENIGN                  --
##---------------------------------------------------------------


ref=DEenrichRPlot(reference,
                  ident.1 = "OVSAHO",
                  max.genes = 100,
                  test.use = "wilcox",
                  p.val.cutoff = 0.05,
                  cols = NULL,
                  enrich.database = "MSigDB_Oncogenic_Signatures",
                  num.pathway = 30,
                  logfc.threshold = 0.1,
                  return.gene.list = T)

query=DEenrichRPlot(ovarian.query, 
                    ident.1 = "OVSAHO",
                    max.genes = 100,
                    test.use = "wilcox",
                    p.val.cutoff = 0.05,
                    cols = NULL,
                    enrich.database = "MSigDB_Oncogenic_Signatures",
                    num.pathway = 30,
                    logfc.threshold = 0.1,
                    return.gene.list = T)



venpos<- venndetail(list(reference =ref$pos$MSigDB_Oncogenic_Signatures.Term , 
                         Query = query$pos$MSigDB_Oncogenic_Signatures.Term
))
venneg<- venndetail(list(reference =ref$neg$MSigDB_Oncogenic_Signatures.Term ,
                         Query = query$neg$MSigDB_Oncogenic_Signatures.Term
))


allpos=getSet(venpos, subset = c("Shared"))

allpos %>% 
  mutate(ID="Pos") %>% 
  mutate(Status="BENIGN") %>% 
  mutate(CellLine="OVSAHO") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->allposOVSAHO

allneg=getSet(venneg, subset = c("Shared"))  
allneg %>% 
  mutate(ID="Neg") %>% 
  mutate(Status="BENIGN") %>% 
  mutate(CellLine="OVSAHO") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->allnegOVSAHO


OVSAHO_BENIGN <- rbind(allposOVSAHO, allnegOVSAHO)

#Extract genes as well for comparison to expression analysis
query$pos %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryposgenes

allposOVSAHO %>% 
  left_join(queryposgenes, by="Pathway")-> OVSAHO_pos_benign

query$neg %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryneggenes

allnegOVSAHO %>% 
  left_join(queryneggenes, by="Pathway")-> OVSAHO_neg_benign

OVSAHO_GENES_BENIGN<- rbind(OVSAHO_neg_benign, OVSAHO_pos_benign)

##----------------------------------------------------------------
##                 OVSAHO Cell Line - NO BENIGN                 --
##----------------------------------------------------------------


ref=DEenrichRPlot(reference,
                  ident.1 = "OVSAHO",
                  max.genes = 100,
                  test.use = "wilcox",
                  p.val.cutoff = 0.05,
                  cols = NULL,
                  enrich.database = "MSigDB_Oncogenic_Signatures",
                  num.pathway = 30,
                  logfc.threshold = 0.1,
                  return.gene.list = T)

query.tumors=DEenrichRPlot(ovarian.tumors, 
                           ident.1 = "OVSAHO",
                           max.genes = 100,
                           test.use = "wilcox",
                           p.val.cutoff = 0.05,
                           cols = NULL,
                           enrich.database = "MSigDB_Oncogenic_Signatures",
                           num.pathway = 30,
                           logfc.threshold = 0.1,
                           return.gene.list = T)

ven.tumors.pos <- venndetail(list(reference =ref$pos$MSigDB_Oncogenic_Signatures.Term , 
                      Query = query.tumors$pos$MSigDB_Oncogenic_Signatures.Term
))
ven.tumors.neg <- venndetail(list(reference =ref$neg$MSigDB_Oncogenic_Signatures.Term , 
                      Query = query.tumors$neg$MSigDB_Oncogenic_Signatures.Term
))

#Do the same for Tumour data
tumors.pos=getSet(ven.tumors.pos, subset = c("Shared"))

tumors.pos %>% 
  mutate(ID="Pos") %>% 
  mutate(Status="NO_BENIGN") %>% 
  mutate(CellLine="OVSAHO") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->tumors.posOVSAHO

tumors.neg=getSet(ven.tumors.neg, subset = c("Shared"))

tumors.neg %>% 
  mutate(ID="Neg") %>% 
  mutate(Status="NO_BENIGN") %>% 
  mutate(CellLine="OVSAHO") %>% 
  select(!Subset) %>% 
  rename(Pathway=Detail)->tumors.negOVSAHO

OVSAHO_NOBENIGN <- rbind(tumors.posOVSAHO, tumors.negOVSAHO)
OVSAHO=rbind(OVSAHO_BENIGN, OVSAHO_NOBENIGN)

#Extract genes as well for comparison to expression analysis
query.tumors$pos %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term, 
           MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryposgenes

tumors.posOVSAHO %>% 
  left_join(queryposgenes, by="Pathway")-> OVSAHO_pos_nobenign

query.tumors$neg %>% 
  select(c(MSigDB_Oncogenic_Signatures.Term,
           MSigDB_Oncogenic_Signatures.Genes)) %>% 
  rename(Pathway=MSigDB_Oncogenic_Signatures.Term) %>% 
  rename(query_Genes=MSigDB_Oncogenic_Signatures.Genes)-> queryneggenes

tumors.negOVSAHO %>% 
  left_join(queryneggenes, by="Pathway")-> OVSAHO_neg_nobenign

OVSAHO_GENES_NOBENIGN<- rbind(OVSAHO_neg_nobenign, OVSAHO_pos_nobenign)
OVSAHO_GENES<-rbind(OVSAHO_GENES_NOBENIGN, OVSAHO_GENES_BENIGN)

#Check which ones have conflicting results
OVSAHO %>%
  group_by(Pathway, CellLine, Status) %>% 
  count() %>% 
  filter(n>1)

OVSAHO_alt <- OVSAHO
OVSAHO_alt %>% 
  pivot_wider(names_from=Status, values_from=ID) %>% 
  mutate(Pathway_bold=ifelse(BENIGN==NO_BENIGN, 
                             paste("**", Pathway, "**", sep=""), Pathway)) %>% 
  mutate(Pathway_bold=ifelse(is.na(Pathway_bold), Pathway, Pathway_bold)) %>% 
  mutate(Pathway=Pathway_bold) %>% 
  select(!Pathway_bold) %>% 
  pivot_longer(c(BENIGN, NO_BENIGN), names_to="Status", values_to="ID") %>% 
  drop_na()-> OVSAHO_plot

OVSAHO_plot

# Plot pathways from OVSAHO for benign and no benign
p3 <- ggplot(OVSAHO_plot, aes(x=Status , y=Pathway, colour = ID)) +
  geom_point(alpha=2,size = 6)+
  theme_classic()+
  theme(axis.text.y = element_markdown(),plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  ggtitle("OVSAHO")

p3

# Common pathways between benign and no benign cells
par(mfrow=c(1,3))
p1+p2+p3

message <- "Common genes."

system2(command = "PowerShell", 
        args = c("-Command", 
                 "\"Add-Type -AssemblyName System.Speech;",
                 "$speak = New-Object System.Speech.Synthesis.SpeechSynthesizer;",
                 paste0("$speak.Speak('", message, "');\"")
        ))

##===============================================================
##  SECTION 5: Find common genes between the shared pathways   ==
##  (benigh and no benign) and all markers from each cluster   ==
##===============================================================

# Merge the different Pathway dataframes for each cell line
pathways<-rbind(KURAM_GENES, CAOV3_GENES, OVSAHO_GENES)

# Split genes by ";"
split_genes <- strsplit(pathways$query_Genes, ";")

# Structured representation of the data for further analysis or visualization
pathways_exp <- data.frame(
  Pathway = rep(pathways$Pathway, sapply(split_genes, length)),
  ID = rep(pathways$ID, sapply(split_genes, length)),
  Status = rep(pathways$Status, sapply(split_genes, length)),
  CellLine = rep(pathways$CellLine, sapply(split_genes, length)),
  Genes = unlist(split_genes)
)

# Import Gene Expression from a csv file
expression <- read.csv("all_markers_proteome.csv")

# Extract top50 genes grouped by cluster
expression %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),group_by=TRUE) %>% 
  top_n(avg_log2FC,n=50) %>% 
  ungroup() -> top50

# Find common genes using inner join by genes and cluster
pathways_exp %>% 
  inner_join(expression, by=c("Genes"="gene", "CellLine"="cluster")) -> shared

# pathways_exp %>% 
#   inner_join(top50, by=c("Genes"="gene", "CellLine"="cluster"))-> shared
# only 5 and all in KURAM

# A summary of shared pathways by cell line
shared %>% 
  group_by(Pathway, CellLine) %>% 
  count() %>% 
  arrange(desc(n))

# Write the shared genes to a csv file
shared %>% 
  write_csv("pathway_exp_commongenes.csv")

# Import table with common genes between this and pathway analysis
common_genes<-read.csv("pathway_exp_commongenes.csv")

# Heatmap for the common genes between common pathways and all markers
P2=DoHeatmap(ovarian.query,features = common_genes$Genes,group.bar = 
               TRUE,label=TRUE,angle=30,size = 3.5)+  guides(colour=FALSE)
P2

##================================================================
##                      END OF THE SCRIPT                       ==
##================================================================


# This is the end of the script
message <- "This is the end of the script. 
            Please contact the author for questions or doubts"

system2(command = "PowerShell", 
        args = c("-Command", 
                 "\"Add-Type -AssemblyName System.Speech;",
                 "$speak = New-Object System.Speech.Synthesis.SpeechSynthesizer;",
                 paste0("$speak.Speak('", message, "');\"")
        ))


