#First attempt at phylosignal analysis
#Calculating K-statistic(Blomberg) and p-value

library(phytools)
library(DESeq2)


#set up

DATAPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-10-K-statistic-calculation/data/"
FIGUREPATH = "~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-02-10-K-statistic-calculation/figures/"
odd <- seq(1,26,2)
even <- seq(2, 26, 2)

#KSNP tree provided by Ali:

full_tree <- read.tree(text ="((HM17:0.00428,HM57:0.00876)1.000:0.08988,(HM65:0.00824,HM68:0.00963)1.000:0.09530,(((HM54:0.08292,HM66:0.07434)1.000:0.04338,(HM43:0.11275,((HM60:0.24376,((HM7:0.07135,(HM14:0.08006,HM3:0.07130)1.000:0.02392)1.000:0.40723,(HM46:0.23512,(HM26:0.19620,(HM1:0.01633,HM69:0.01337)1.000:0.17697)1.000:0.04329)1.000:0.10117)1.000:0.07596)1.000:0.33315,(HM56:0.08289,HM6:0.08508)1.000:0.08846)1.000:0.04832)1.000:0.01448)1.000:0.01347,(HM86:0.10068,HM27:0.11440)0.970:0.01392)1.000:0.01738);")
#focusing on 13 strains 
tree <- drop.tip(full_tree, c("HM65", "HM27", "HM26", "HM69", "HM46", "HM60"))

#working with normalized counts for the core
#NOTE: might want to do this with transformed values?
#NOTE: some issues with these counts, i.e. paralogs
norm_counts <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2017-01-30-htseq-counts/normalized_counts_core.csv", header = TRUE, row.names = 1)

#URINE vs UTI
norm_urine_counts <- norm_counts[,odd]
norm_uti_counts <- norm_counts[,even]


#calculating K statistic, saving values of K-statistic P-value and normalized counts
calculateKForDataFrame <- function(trait_matrix, tree, filename) {
        
        K_values <- c()
        P_values <- c()
        genes_with_signal <- c()
        mat <- matrix(NA, nrow=nrow(trait_matrix), ncol=3)
        
        for (i in 1:nrow(trait_matrix)){
     
                gene_name = rownames(trait_matrix[i,])
                expression_vector <- as.numeric(trait_matrix[i,])
                names(expression_vector) <- c("HM1", "HM3", "HM6", "HM7", "HM14", "HM17", "HM43", "HM54", "HM56", "HM57", "HM66", "HM68", "HM86")
                K <- phylosig(tree, expression_vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                         control=list())$K
                P <- phylosig(tree, expression_vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                              control=list())$P
                if (P <= 0.05){
                        genes_with_signal <- append(genes_with_signal, c(gene_name, K, P))
                        mat[i,]<- c(gene_name, K, P)
                }
                K_values <- append(K_values, K)
                P_values <- append(P_values, P)
        }
        
        mat <- mat[complete.cases(mat),]
       
        phylo_sig <- mat[,c(2,3)]
        rownames(phylo_sig) <- mat[,1]
        signal_names <- mat[,1]
        colnames(phylo_sig)<- c("K-statistic", "P-value")
        phylosig_traits <- trait_matrix[signal_names,]
        
        phy <- list("phylosig_K" = as.data.frame(phylo_sig), "phylosig_gene_counts" = phylosig_traits)
        fn = paste0(DATAPATH, filename)
        write.csv(phylo_sig, fn, row.names=TRUE)
        return (phy)
}


urine <- calculateKForDataFrame(norm_urine_counts, tree, "core_genes_urine_K_statistic.csv")
uti <- calculateKForDataFrame(norm_uti_counts, tree, "core_genes_uti_K_statistic.csv")


core_genes_with_physignal <- merge(urine$phylosig_gene_counts, uti$phylosig_gene_counts, by="row.names", all=TRUE)
write.csv(core_genes_with_physignal, paste0(DATAPATH, "core_genes_with_physignal.csv"), row.names = TRUE)

#look at the intersection of urine and uti datasets

strong_signal <- core_genes_with_physignal[complete.cases(core_genes_with_physignal),]




names <- strong_signal[,1]# should be 51 names

mat <- assay(core_dds_transformed)[names,]
head (mat)
mat <- mat-rowMeans(mat)
df <- as.data.frame(colData(core_dds_transformed)[,c("MEDIA", "PRED_PHYLO")])
pheatmap(mat, annotation_col=df)


#URINE
ur_p <- urine_phylosig[order(urine_phylosig$K.statistic, decreasing = TRUE),]
ur_p <- ur_p[ur_p$P.value <= 0.001,]

ur_names <- rownames(ur_p)

mat_ur <- assay(core_urine_transformed)[ur_names,]
mat_ur <- mat_ur-rowMeans(mat_ur)
df <- as.data.frame(colData(core_urine_transformed)[,c("PRED_PHYLO", "STRAIN")])
pheatmap(mat_ur, annotation_col=df)
