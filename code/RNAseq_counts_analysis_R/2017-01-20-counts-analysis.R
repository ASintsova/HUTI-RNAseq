#Dependencies
library(DESeq2)
library (pheatmap)
library(RColorBrewer)
library (ggplot2)
source("~/git_repos/HUTI-RNAseq/code/RNAseq_counts_analysis_R/RNAseq_analysis_utilities.R")
#Importing the counts
combined_counts <- read.csv("~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/combined_counts.csv", header = TRUE, row.names = 1)

#Importing metadata
experiment_info <- read.csv("~/git_repos/HUTI-RNAseq/analysis/Patient_meta_info.csv", header = TRUE, row.names =1)

head(combined_counts)

#Isolating Core Genes

core_genes = combined_counts[complete.cases(combined_counts),]
summary(core_genes)
(even <- seq(2, 26, 2))
(odd <- seq(1, 26, 2))

#Seperationg urine and uti samples
core_urine <- core_genes[,odd]
core_uti <- core_genes[,even]
urine_info <- experiment_info[odd, ]
uti_info <- experiment_info[even,]
write.csv(core_genes, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/core_gene_counts.csv", row.names=TRUE)
write.csv(core_urine, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/core_counts_urine.csv", row.names=TRUE)
write.csv(core_uti, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/core_counts_uti.csv", row.names=TRUE)
write.csv(urine_info, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/urine_info.csv", row.names=TRUE)
write.csv(uti_info, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/uti_info.csv", row.names=TRUE)



# Construct DESeqDataSet object from the matrix of counts and the experiment_info
(core_dds <- DESeqDataSetFromMatrix(countData = as.matrix(core_genes), colData = experiment_info,
                                     design = ~MEDIA))
(core_uti_dds <- DESeqDataSetFromMatrix(countData = as.matrix(core_uti), colData = uti_info,
                                        design = ~PRED_PHYLO))
(core_urine_dds <- DESeqDataSetFromMatrix(countData = as.matrix(core_urine), colData = urine_info,
                                        design = ~PRED_PHYLO))
#Pre-filtering the dataset: removing rows of DESeqDataSet that have not counts
#or only single count across all samples


core_dds_filtered <- core_dds[rowSums(counts(core_dds)) > 1,]
dim(core_dds_filtered)

core_urine_dds_filtered <- core_urine_dds[rowSums(counts(core_urine_dds))>1, ]
dim(core_urine_dds_filtered)
core_uti_dds_filtered <- core_uti_dds[rowSums(counts(core_uti_dds)) >1,]
dim(core_uti_dds_filtered)

## COMMENT: This filtering step did not help, probably want to filter only based on UTI samples, UR samples will
## probalby have counts, or need to set a differnt threshhold?  

#Data transform: DESeq2 specific, there are disadvantages to doing clustering analysis on RPKMs,
# and even normal log transforms
core_dds_transformed <- rlog(core_dds, blind=FALSE)#FALSE settin gwill means differences due to 
# design are not taken into account

par(mfrow = c(2, 7))

for (i in odd){
        print (i)
        print (i+1)
        plot(assay(core_dds_transformed)[,i:(i+1)], pch =16, cex =0.3)
}
# How do I save this?


#R function dist to calculate Euclidean distance between samples
#Transform matrix
core_dists <- dist(t(assay(core_dds_transformed)))
dists_matrix <- as.matrix(core_dists)
rownames(dists_matrix)<- paste(core_dds_transformed$STRAIN, core_dds_transformed$MEDIA, sep = "_")
colnames(dists_matrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(dists_matrix,
         clustering_distance_rows=core_dists,
         clustering_distance_cols=core_dists,
         col=colors)


# Look at other clustering methods: hclust??

plotPCA(core_dds_transformed, intgroup = c("MEDIA"))

#for (name in colnames(experiment_info)){
#        plotPCA(core_dds_transformed, intgroup = c(name))
#}
# this did not work, need to see how I can generate and save them somewhere else


#PCA using ggplot

pcaData <- plotPCA(core_dds_transformed, intgroup = c("MEDIA", "ANCESTRY"), returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

ggplot(data = pcaData, aes(x=PC1, y=PC2, color = MEDIA, shape = ANCESTRY)) + 
        geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% var")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed()

core_DE <- DESeq(core_dds)
#Set LFC threshold 1 
core_DE_results <- results(core_DE, lfcThreshold = 1)
core_DE_results_sign <- subset(core_DE_results, padj < 0.05)
head(core_DE_results_sign[order(core_DE_results_sign$log2FoldChange, decreasing = TRUE),], 20)

core_DE_genes_final <- core_DE_results_sign[order(core_DE_results_sign$log2FoldChange, decreasing = TRUE),]

topGene <- rownames(core_DE_results_sign)[which.min(core_DE_results_sign$padj)]
plotCounts(core_dds, gene = topGene, intgroup = c("MEDIA"))

for (i in 1:14){
      
        print (rownames(top)[i])
        plotCounts(core_dds, gene = rownames(top)[i], intgroup = c("MEDIA"))
}

gene_names <- rownames(core_DE_genes_final)

gene_names_edited <- list()
for (name in gene_names){
        print (name)
        new_name <- unlist(strsplit(name, "_"))[2
        print(new_name)
        gene_names_edited <- append(new_name, gene_names_edited)
}

gene_names_edited <- unlist(gene_names_edited)

