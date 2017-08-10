#Brainstroming analysis of transcription regulatory network in my E.coli data


# First want to look at whether there are differenctially expressed genes between different genotypes

data <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/data/core_gene_counts.csv",
                 row.names =1)
syms <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/data/DE-seq-results-all-core-edited.csv",
                 row.names =1)[, 7]

data2 <- cbind(data, syms)
info <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/data/Patient_meta_info.csv",
                 row.names =1)

#urine counts
odd <-seq(1,26,2)
even <- seq(2,26,2)
data3 <- data2[,even]
info3 <- info[even,]


#remove E phylogroup from analysis (sample UR1)
data4 <- data3[,2:13]
info4 <- info3[2:13,]


# now will perform DE based on phylogroup

library(DESeq2)
(dds <- DESeqDataSetFromMatrix(countData = as.matrix(data4), colData = info4,
                                    design = ~PRED_PHYLO))

ddsT <- rlog(dds, blind=FALSE)

plotPCA(ddsT, intgroup = c("PRED_PHYLO"))


DE <- DESeq(dds)

de_results <- results(DE, lfcThreshold = 0.5)
                    
de_sign <- subset(de_results, padj < 0.1)

de_sign <- de_sign[order(de_sign$log2FoldChange, decreasing = TRUE),]


names <- rownames(de_sign)


#get gene symbols for these

data5 <- data2[rownames(data2) %in% names,]
rownames(de_sign) <- data5[,27]
de_names <- data5[,27]




urine <- cbind(as.data.frame(de_sign), de_names)
#de_names are all messed up

write.csv(urine[,1:6], "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-07-25/data/2017-07-25-urine-by-phyG.csv")
uti <- cbind(as.data.frame(de_sign), de_names) 
write.csv(uti[,1:6], "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-07-25/data/2017-07-25-uti-by-phyG.csv")

both_names <- intersect(rownames(urine), rownames(uti))

both <- cbind(urine[rownames(urine) %in% both_names,1:6], uti[rownames(uti) %in% both_names,1:6])

dim(de_sign)
length(data5[,27])




# now want to repeat the same for UTI
#go up change odd to even


