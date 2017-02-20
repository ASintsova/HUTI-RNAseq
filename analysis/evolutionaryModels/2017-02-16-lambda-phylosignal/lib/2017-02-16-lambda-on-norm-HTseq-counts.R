#Calculating lambda on the same dataset

library(phytools)
library(DESeq2)
library(pheatmap)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/KStatFromTraitMatrix.R")
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LambdaFromTraitMatrix.R")
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")


DATAPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-16-lambda-phylosignal/data/"
FIGUREPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-16-lambda-phylosignal/figures/"


all_data <- LoadHTSeqCountsandTree()
plotTree(all_data$tree)

#run lambda calculations on log2 transformed urine and uti counts
lambda_urine <- calculateLambdaForDataFrame(all_data$urine_counts_l, all_data$tree,
                                            "2017-02-16-lambda-on-urine.csv")
#this found 240 genes
lambda_uti <- calculateLambdaForDataFrame(all_data$uti_counts_l, all_data$tree,
                                            "2017-02-16-lambda-on-uti.csv")
#found 70 genes, consistent with K stat calculations

#combine to find intersection

lambda_core <- merge(lambda_urine$phylosig_gene_counts, lambda_uti$phylosig_gene_counts, 
                     by ="row.names", all = TRUE )


strong_signal <- lambda_core[complete.cases(lambda_core),]
row.names(strong_signal) <- strong_signal[,1]
strong_signal <- strong_signal[,2:27]
dim(strong_signal)
annot <- all_data$exp_info[, c("MEDIA", "PRED_PHYLO")]

pheatmap(strong_signal, annotation_col = annot, cellwidth = 20, cellheight = 15,
         #display_numbers = counts, cutree_cols = 2, cutree_rows = 4,
         filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-16-lambda-phylosignal/figures/2017-02-16-core-phylosig-lambda.png")

#urine
urine_low_pval <- lambda_urine$phylosig_lambda[lambda_urine$phylosig_lambda$`P-value` == 0.001,]
urine_low_pval_names <- row.names(urine_low_pval)
urine_counts <- round(all_data$urine_counts[urine_low_pval_names,]/1000, 1)
annot <- all_data$urine_info[, c("HISTORY", "PRED_PHYLO")]
pheatmap(all_data$urine_counts_l[urine_low_pval_names,], annotation_col = annot, cellwidth = 20, cellheight = 15,
         display_numbers = urine_counts, cutree_cols = 2, cutree_rows = 4,
         filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-16-lambda-phylosignal/figures/2017-02-16-urine-low-pval-phylosig-lambda.png")




uti_low_pval <- lambda_uti$phylosig_lambda[lambda_uti$phylosig_lambda$`P-value` == 0.001,]
dim(uti_low_pval)
uti_low_pval_names <- row.names(uti_low_pval)
uti_counts <- round(all_data$uti_counts[uti_low_pval_names,]/1000, 1)
annot <- all_data$uti_info[, c("HISTORY", "PRED_PHYLO")]
pheatmap(all_data$uti_counts_l[uti_low_pval_names,], annotation_col = annot, cellwidth = 20, cellheight = 15,
         display_numbers = uti_counts, cutree_cols = 2, cutree_rows = 4,
         filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-16-lambda-phylosignal/figures/2017-02-16-uti-low-pval-phylosig-lambda.png")


pheatmap(lambda_uti$phylosig_gene_counts, annotation_col = annot, cellwidth = 15, cellheight = 8, cutree_cols = 2, 
         filefilename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-16-lambda-phylosignal/figures/2017-02-16-uti-all-phylosig-lambda.png")

