#source("https://bioconductor.org/biocLite.R")
library(DESeq2)

results = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/differential_expression_analysis/"
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
process_dds <- function (dds, prefix){
 
  analysis <- DESeq(dds)
  res <- results(analysis)
  resSig <- subset(res, (padj < 0.05)& (log2FoldChange > 2 | log2FoldChange < -2))
  write.csv(resSig[order(resSig$log2FoldChange, decreasing = TRUE),], paste0(results, prefix, "_significant.csv"),
            row.names = TRUE, quote = FALSE)
  write.csv(res[order(res$log2FoldChange, decreasing = TRUE),], paste0(results, prefix, "_all.csv"),
            row.names = TRUE, quote = FALSE)
  esf <- estimateSizeFactors(analysis)
  norm_counts <- counts(esf, normalized=TRUE)
  write.csv(norm_counts, paste0(results, prefix, "_DSeq2_norm_counts.csv"), quote = FALSE)
  output <-  list("results" = resSig[order(resSig$log2FoldChange, decreasing = TRUE),], "dds"= analysis)
  return (output)
}

# Differential Expression Analysis by 'MEDIA'

dds_media <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = experiment_info, design = ~MEDIA)
dds_media$MEDIA <- relevel(dds_media$MEDIA, ref = "URINE")
media <- process_dds(dds_media, "media_de_genes")

# Differential Expression Analysis by phylogroup, controlling for 'MEDIA'
# Need to drop E
drop <- colnames(counts) %in% c("HM01_UR", "HM01_UTI")
counts_pg <- counts[!drop]
experiment_info_pg <- experiment_info[!drop, ]

dds_pg <- DESeqDataSetFromMatrix(countData = as.matrix(counts_pg), colData = experiment_info_pg, 
                              design = ~MEDIA + PRED_PHYLO)
pg_controlled <- process_dds(dds_pg, "pg_contr_de_genes" )


# Drop all of the UTI
drop_uti = grepl("UTI", colnames(counts_pg))
counts_urine = counts_pg[!drop_uti]
experiment_info_urine <- experiment_info_pg[!drop_uti, ]

dds_ur <- DESeqDataSetFromMatrix(countData = as.matrix(counts_urine), colData = experiment_info_urine, 
                              design = ~PRED_PHYLO)
pg_urine <- process_dds(dds_ur, "pg_urine_de_genes")


# Drop all of the UR samples
drop_ur = grepl("UR", colnames(counts_pg))
counts_uti = counts_pg[!drop_ur]
experiment_info_uti <- experiment_info_pg[!drop_ur, ]

dds_uti <- DESeqDataSetFromMatrix(countData = as.matrix(counts_uti), colData = experiment_info_uti, 
                                 design = ~PRED_PHYLO)
pg_uti <- process_dds(dds_uti, "pg_uti_de_genes")