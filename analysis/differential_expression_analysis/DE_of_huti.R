source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library(DESeq2)
library (pheatmap)
library(RColorBrewer)
library (ggplot2)
library(dplyr)

counts_file <- "/Users/annasintsova/git_repos/HUTI-RNAseq/results/methods/core_counts_and_tpms.csv"
colData_file <- "/Users/annasintsova/git_repos/HUTI-RNAseq/data/huti_patient_info_for_DE.csv"

counts <- read.table(counts_file, row.names =1, sep = ",", header = TRUE)

old_names <- c()
new_names <- c()
for (c in colnames(counts)){
  if (grepl( "count", c)){
    old_names <-append(old_names, c)
    nn <- unlist(strsplit(c, "_count")[[1]])
    new_names <- c(new_names, nn)
  }
}
counts <- counts[old_names]
colnames(counts) <- new_names

experiment_info <- read.csv(colData_file, header = TRUE, row.names =1)


# Differential Expression Analysis by 'MEDIA'

dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = experiment_info, design = ~MEDIA)

