##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Date: Mon Jun 5 10:28:18 2024                                                                    ::
##  File Name: Stats.R                                                                               ::
##  Author: Eleni Theofania Skorda                                                                   ::
##                                                                                                   ::
##  Description:                                                                                     ::
##      This script was was used to calculate differential expression analysis statistics.           ::
##      Then we performed differential expression analysis and lastly we run the OmicLoupe           ::
##      through Shiny App.                                                                           ::
##                                                                                                   ::
##  List of Functions:                                                                               ::
##      none                                                                                         ::
##                                                                                                   ::
##  Procedure:                                                                                       ::
##      1. Prepare datasets to a proper format                                                       ::
##      2. Run NormalyzerDE without Log Transformation because the dataset is already transformed.   ::
##      3. Perform differential expression analysis                                                  ::
##      4. Run OmicLoupe through ShinyApp                                                            ::
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Load necessary libraries
library(NormalyzerDE)
library(OmnicLoupe)
library(ggplot2)
library(ggrepel)

# Prepare datasets to a proper format 
peptide <- read.delim("ratio_peptide_MD.tsv", sep = "\t",check.names = FALSE)

exprBT = peptide[,-c(1,3:6)]
write.table(exprBT, file = "peptide_matrix.tsv", sep = "\t", quote = FALSE, 
            row.names = FALSE)


## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE, dev="png", message=FALSE,
                      error=FALSE, warning=TRUE)
## -----------------------------------------------------------------------------

# Run NormalyzerDE
library(NormalyzerDE)
normalyzerDE(jobName="peptide", designPath="table_design.tsv", 
             dataPath="peptide_matrix.tsv", comparisons=c("Tumor-Normal", 
                                                          "Normal-Not_reported",
                                                          "Tumor-Not_reported"))

# Perform differential expression analysis
stats <- read.delim("peptide_stats.tsv", sep = "\t",check.names = FALSE)
de <- tmp[complete.cases(stats), ]

# Rename columns
colnames(de)[2] ="pvalue"
colnames(de)[8] ="log2FoldChange"

# Plot de genes
ggplot(data=de, aes(x=log2FoldChange, y=pvalue)) + geom_point()

# Convert directly in the aes()
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
p
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + 
  theme_minimal()
p
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

# Add a column of NAs
de$diffexpressed <- "NO"


# Set limits
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.6 & de$pvalue < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.6 & de$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed))
+ geom_point() + theme_minimal()

# Add lines as before
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# Automation

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Create new column which will contain differential expressed genes
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$Gene[de$diffexpressed != "NO"]

#Plot it
ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed,
                    label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# Organize the labels in a more nice way
ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, 
                    label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# Run OmicLoupe
OmicLoupe::runApp()
