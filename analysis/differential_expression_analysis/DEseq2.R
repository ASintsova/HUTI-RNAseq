library(DESeq2)
library(topGO)
library("org.EcK12.eg.db")

cts <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/core_counts.csv", header = TRUE, row.names = 1)

coldata <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/coldata.csv", header = TRUE, row.names = 1)

#columns of cts need to be in the same order as rows of coldata!!!
all(rownames(coldata) == colnames(cts)) # Needs to be True



# Analysis using all samples

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment)

#Differential Expression Analysis

analysis <- DESeq(dds)
res <- results(analysis)
res05 <- results(analysis, alpha=0.05)
res05_ordered <- res05[order(res05$log2FoldChange),]
write.csv(res05_ordered, "~/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/all_strains_DEseq.csv", row.names=TRUE)

plotMA(res05, ylim=c(-2,2))

#Using only 'the best'
# Best: ["HM56", "HM14", "HM43", "HM54", "HM86"]

best <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/best_counts.csv", header = TRUE, row.names = 1)
coldata_best <- coldata[colnames(best),]
dds_best <- DESeqDataSetFromMatrix(countData = best,
                              colData = coldata_best,
                              design = ~ treatment)

analysis_best <- DESeq(dds_best)
res05_best <- results(analysis_best, alpha=0.05)
res05_best_ordered <- res05_best[order(res05_best$log2FoldChange),]
write.csv(res05_best_ordered, "~/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/best_strains_DEseq.csv", row.names=TRUE)


resLFC <- lfcShrink(analysis_best, coef=2)
plotM(resLFC, ylim=c(-3,3))





source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("graph")
biocLite("Rgraphviz")
biocLite("SparseM")
biocLite("GOFunction")


library("GOFunction")
data(exampledata)

sigTerm <- GOFunction(interestGenes, refGenes, organism="org.Hs.eg.db",
                      ontology="BP", fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05,
                      poth=0.05, peth=0.05, bmpSize=2000, filename="sigTerm")



ig_df <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/best_gene_of_interst_eg.csv", header = TRUE, row.names = 1)
ig <- rownames(ig_df)

ref_df <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/best_refGene_eg.csv", header = TRUE, row.names = 1)
rg <- rownames(ref_df)


