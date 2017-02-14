# Second attempt, going to look at UTI samples alone
# I want to look of results differe if I use normalized counts vs. log transformed counts,
# also look at different visualization options

# set up at of 2017-02-10
# also import a function calculateKForDataFrame that calculates K statistic given a trait matrix, a tree,
# a filename for the output

library(phytools)
library(DESeq2)
library(pheatmap)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/KStatFromTraitMatrix.R")


DATAPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-14-K-statistic-contd/data/"
FIGUREPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-10-K-statistic-contd/figures/"
odd <- seq(1,26,2)
even <- seq(2, 26, 2)
#metadata
uti_info <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-01-14-htseq-counts/data/uti_info.csv", header = TRUE, row.names=1)



#KSNP tree provided by Ali:

full_tree <- read.tree(text ="((HM17:0.00428,HM57:0.00876)1.000:0.08988,(HM65:0.00824,HM68:0.00963)1.000:0.09530,(((HM54:0.08292,HM66:0.07434)1.000:0.04338,(HM43:0.11275,((HM60:0.24376,((HM7:0.07135,(HM14:0.08006,HM3:0.07130)1.000:0.02392)1.000:0.40723,(HM46:0.23512,(HM26:0.19620,(HM1:0.01633,HM69:0.01337)1.000:0.17697)1.000:0.04329)1.000:0.10117)1.000:0.07596)1.000:0.33315,(HM56:0.08289,HM6:0.08508)1.000:0.08846)1.000:0.04832)1.000:0.01448)1.000:0.01347,(HM86:0.10068,HM27:0.11440)0.970:0.01392)1.000:0.01738);")
#focusing on 13 strains 
tree <- drop.tip(full_tree, c("HM65", "HM27", "HM26", "HM69", "HM46", "HM60"))

#import normalized counts

norm_counts <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2017-01-30-htseq-counts/normalized_counts_core.csv", header = TRUE, row.names = 1)

#UTI
norm_uti_counts <- norm_counts[,even]

#transform

uti_log <- log2(norm_uti_counts + 1)

signal_on_counts <- calculateKForDataFrame(norm_uti_counts, tree, "uti_signal_on_counts.csv")
signal_on_log_transform <- calculateKForDataFrame(uti_log, tree, "uti_signal_on_log_transform.csv")
# pretty similar output

#order each by P-value

count_phylosig <- signal_on_counts$phylosig_K
log_phylosig <- signal_on_log_transform$phylosig_K

dim(count_phylosig)
dim(log_phylosig)

sig_count_phylosig<- count_phylosig[order(count_phylosig$`K-statistic`, decreasing = TRUE),]
sig_log_phylosig<- log_phylosig[order(log_phylosig$`K-statistic`, decreasing = TRUE),]

uti <- signal_on_counts$phylosig_gene_counts
uti_mean_norm <- uti - rowMeans(uti)
df <- uti_info[, c("PRED_PHYLO", "HISTORY")]
pheatmap(uti, annotation_col=df) #this doesn't look very good because spread is too big
#does not matter whether its centered on row mean

uti_log2 <- signal_on_log_transform$phylosig_gene_counts
uti_log2_mean_norm <- uti_log2 - rowMeans(uti_log)
rownames(df) == colnames(uti_log2)
pheatmap(uti_log2, annotation_col = df) #looks ok, clusters based on phylogroup
pheatmap(uti_log2_mean_norm, annotation_col=df)# does not look quite as good, a bit of a mess
#mat_ur <- assay(core_urine_transformed)[ur_names,]
#mat_ur <- mat_ur-rowMeans(mat_ur)
#df <- as.data.frame(colData(core_urine_transformed)[,c("PRED_PHYLO", "STRAIN")])
pheatmap(mat_ur, annotation_col=df)


