#Repeating yesterdays analysis for UR and combined data
#Hope to improve visualization for meeting for tomoroow

library(phytools)
library(DESeq2)
library(pheatmap)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/KStatFromTraitMatrix.R")


DATAPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-15-K-statistic-UR-and-combined/data/"
FIGUREPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-15-K-statistic-UR-and-combined/figures/"

odd <- seq(1,26,2)
even <- seq(2, 26, 2)
#metadata
urine_info <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-01-14-htseq-counts/data/urine_info.csv", header = TRUE, row.names=1)
experiment_info <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/data/Patient_meta_info.csv", header = TRUE, row.names=1)


#KSNP tree provided by Ali:

full_tree <- read.tree(text ="((HM17:0.00428,HM57:0.00876)1.000:0.08988,(HM65:0.00824,HM68:0.00963)1.000:0.09530,(((HM54:0.08292,HM66:0.07434)1.000:0.04338,(HM43:0.11275,((HM60:0.24376,((HM7:0.07135,(HM14:0.08006,HM3:0.07130)1.000:0.02392)1.000:0.40723,(HM46:0.23512,(HM26:0.19620,(HM1:0.01633,HM69:0.01337)1.000:0.17697)1.000:0.04329)1.000:0.10117)1.000:0.07596)1.000:0.33315,(HM56:0.08289,HM6:0.08508)1.000:0.08846)1.000:0.04832)1.000:0.01448)1.000:0.01347,(HM86:0.10068,HM27:0.11440)0.970:0.01392)1.000:0.01738);")
#focusing on 13 strains 
tree <- drop.tip(full_tree, c("HM65", "HM27", "HM26", "HM69", "HM46", "HM60"))


#working with normalized counts for the core
norm_counts <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2017-01-30-htseq-counts/normalized_counts_core.csv", header = TRUE, row.names = 1)
urine_counts <- norm_counts[, odd]
uti_counts <- norm_counts[,even]

#normalized all counts and urine counts
norm_counts_l <- log2(norm_counts +1)
urine_counts_l <- log2(urine_counts +1)
uti_counts_l <- log2(uti_counts +1)
signal_on_uti_counts <- calculateKForDataFrame(uti_counts, tree, "2017-02-15-signal-on-uti-counts.csv")
signal_on_urine_counts <- calculateKForDataFrame(urine_counts, tree, "2017-02-15-signal-on-urine-counts.csv")



#combined ur and uti data
signal_on_core <- merge(signal_on_urine_counts$phylosig_gene_counts, signal_on_uti_counts$phylosig_gene_counts, by="row.names", all=TRUE)

strong_signal <- signal_on_core[complete.cases(signal_on_core),] ## only 15 not 51 since I lowered the P value cut off
rownames(strong_signal) <- strong_signal[,1]
strong_signal <- strong_signal[,2:27]
df <- experiment_info[, c("PRED_PHYLO", "MEDIA")]

log_names <- rownames(strong_signal)
counts <- round(norm_counts[log_names,]/1000,1)
counts <- counts[,c(odd, even)]
dim(strong_signal)
dim(counts)
pheatmap(strong_signal, annotation_col = df, cellwidth = 20, cellheight = 15, 
         display_numbers = counts, cutree_cols = 2, cutree_rows = 4,
         filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-15-K-statistic-UR-and-combined/figures/core_phylosig_low_pval.png")

# urine
urine <- signal_on_urine_counts$phylosig_K
dim(urine)
urine_lp <- urine[as.numeric(as.character(urine$`P-value`)) == 0.001,]
urine <- urine[order(urine$`K-statistic`, decreasing = TRUE),]
ur_names <- rownames(urine)
ur_lp_names <- rownames(urine_lp)

#all urine genes
urine_logs <- urine_counts_l[ur_names,]
dim(urine_logs)
urine_c <- round(urine_counts[ur_names,]/1000,1)
ur_annot <- experiment_info[odd, c("PRED_PHYLO", "HISTORY")]
pheatmap (urine_logs, annotation_col = ur_annot,
          cellwidth = 20, cellheight = 8, display_numbers = urine_c,
          cutree_rows = 3, cutree_cols = 2,
          filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-15-K-statistic-UR-and-combined/figures/urine_phylosig_all.png")

#low pval
urine_logs <- urine_counts_l[ur_lp_names,]
urine_c <- round(urine_counts[ur_lp_names,]/1000,1)
pheatmap (urine_logs, annotation_col = ur_annot,
          cellwidth = 20, cellheight = 8, display_numbers = urine_c,
          cutree_rows = 3, cutree_cols = 2,
          filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-15-K-statistic-UR-and-combined/figures/urine_phylosig_low_pval.png")




#uti

uti <- signal_on_uti_counts$phylosig_K
uti_names <- rownames(uti)


uti_lp <- uti[as.numeric(as.character(uti$`P-value`)) == 0.001,]
uti <- uti[order(uti$`K-statistic`, decreasing = TRUE),]
uti_names_lp <- rownames(uti_lp)

#all uti

uti_logs <- uti_counts_l[uti_names,]
uti_c <- round(uti_counts[uti_names,]/1000,1)
uti_annot <- experiment_info[even, c("PRED_PHYLO", "HISTORY")]
pheatmap (uti_logs, annotation_col = uti_annot,
          cellwidth = 20, cellheight = 20, display_numbers = uti_c,
          cutree_rows = 4, cutree_cols = 2,
          filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-15-K-statistic-UR-and-combined/figures/uti_phylosig_all.png")
#low lp
uti_logs <- uti_counts_l[uti_names_lp,]
uti_c <- round(uti_counts[uti_names_lp,]/1000,1)
uti_annot <- experiment_info[even, c("PRED_PHYLO", "HISTORY")]
pheatmap (uti_logs, annotation_col = uti_annot,
          cellwidth = 20, cellheight = 20, display_numbers = uti_c,
          cutree_rows = 4, cutree_cols = 2,
          filename = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-15-K-statistic-UR-and-combined/figures/uti_phylosig_low_pval.png")





