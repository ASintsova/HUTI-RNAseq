#Calculating lambda on the same dataset

library(phytools)
library(DESeq2)
library(pheatmap)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")


DATAPATH = "~/git_repos/HUTI-RNAseq/analysis/DE/2017-02-17-ribosome-gene-expression/data/"
FIGUREPATH = "~/git_repos/HUTI-RNAseq/analysis/DE/2017-02-17-ribosome-gene-expression/figures/"


all_data <- LoadHTSeqCountsandTree()
ribo_genes <- as.vector(read.table(paste0(DATAPATH, "ribosomal_gene_ids.txt"))[,1])

counts <- all_data$counts
ribo_expr <- counts[ribo_genes,]

sum(is.na(ribo_expr))

ribo_expr <- ribo_expr[complete.cases(ribo_expr),]

ribo_log <- log2(ribo_expr +1)
annot <- all_data$exp_info[, c("MEDIA", "PRED_PHYLO")]
ribo_log_mod <- ribo_log[!rownames(ribo_log) %in% c("15746_rplS", "15749_rpsP"),]#, "17881_ykgO"),]#, "18981_sra", "16752_rpmG",
                                                   # "17882_ykgM", "18114_ribosome-associated_protein" ),]
names <- row.names((ribo_log_mod))
disp_counts <- round(counts[names,]/1000, 1)

pheatmap(ribo_log_mod,annotation_col = annot, display_numbers = disp_counts, 
         cutree_cols = 2, cutree_rows = 4, cellwidth = 20, cellheight = 15, 
         filename = paste0(FIGUREPATH, "2017-02-17-ribosomal-gene-expression.png"))


cm <- colMeans(ribo_log_mod[,seq(1,26,2)])
names(cm)<- c("HM1", "HM3", "HM6", "HM7", "HM14", "HM17", "HM43", "HM54", "HM56", "HM57", "HM66", "HM68", "HM86")
obj <- contMap(all_data$tree, cm)

cm2 <- colMeans(ribo_expr[,seq(1,26,2)])
names(cm2)<- c("HM1", "HM3", "HM6", "HM7", "HM14", "HM17", "HM43", "HM54", "HM56", "HM57", "HM66", "HM68", "HM86")
obj <- contMap(all_data$tree, cm2)

urine_sums <- colSums(ribo_log_mod[,seq(1,26,2)])

mat <- matrix(NA, nrow=nrow(ribo_log_mod), ncol=ncol(ribo_log_mod[,seq(1,26,2)]))
colnames(mat)<- colnames(ribo_expr[,seq(1,26,2)])
rownames(mat) <- row.names(ribo_log_mod)
for (i in 1:nrow(ribo_log_mod)){
       
       print (ribo_log_mod[i, seq(1,26,2)])
        mat[i,] <- as.numeric(round(((ribo_log_mod[i, seq(1,26,2)])/(urine_sums))*100,2))
}




uti_sums <- colSums(ribo_log_mod[,seq(2,26,2)])

mat2 <- matrix(NA, nrow=nrow(ribo_log_mod), ncol=ncol(ribo_log_mod[,seq(2,26,2)]))
colnames(mat2)<- colnames(ribo_expr[,seq(2,26,2)])
rownames(mat2) <- row.names(ribo_log_mod)
for (i in 1:nrow(ribo_log_mod)){
        
        print (ribo_log_mod[i, seq(2,26,2)])
        mat2[i,] <- as.numeric(round(((ribo_log_mod[i, seq(2,26,2)])/(uti_sums))*100,2))
}

ribo_avs_ur <- rowMeans(mat)
ribo_avs_uti <- rowMeans(mat2)
ribos_avs <- matrix(c(ribo_avs_ur, ribo_avs_uti), ncol = 2)
rownames(ribos_avs) <- names(ribo_avs_ur)
ribo_proportions <- merge(mat, mat2, by = "row.names", all =TRUE)
row.names(ribo_proportions) <- ribo_proportions[,1]
ribo_proportions <- ribo_proportions[2:27]
write.csv(ribo_proportions, paste0(DATAPATH, "2017-02-17-proportions-of-ribosomal-gene-expression.csv", row.names = TRUE))



#try the same with normalized counts not with log transforms


urine_sums <- colSums(ribo_expr[,seq(1,26,2)])

mat <- matrix(NA, nrow=nrow(ribo_expr), ncol=ncol(ribo_expr[,seq(1,26,2)]))
colnames(mat)<- colnames(ribo_expr[,seq(1,26,2)])
rownames(mat) <- row.names(ribo_expr)
for (i in 1:nrow(ribo_expr)){
        
        print (ribo_expr[i, seq(1,26,2)])
        mat[i,] <- as.numeric(round(((ribo_expr[i, seq(1,26,2)])/(urine_sums))*100,2))
}




uti_sums <- colSums(ribo_expr[,seq(2,26,2)])

mat2 <- matrix(NA, nrow=nrow(ribo_expr), ncol=ncol(ribo_expr[,seq(2,26,2)]))
colnames(mat2)<- colnames(ribo_expr[,seq(2,26,2)])
rownames(mat2) <- row.names(ribo_expr)
for (i in 1:nrow(ribo_expr)){
        
        print (ribo_expr[i, seq(2,26,2)])
        mat2[i,] <- as.numeric(round(((ribo_expr[i, seq(2,26,2)])/(uti_sums))*100,2))
}

#ribo_avs_ur <- rowMeans(mat)
#ribo_avs_uti <- rowMeans(mat2)
#ribos_avs <- matrix(c(ribo_avs_ur, ribo_avs_uti), ncol = 2)
#rownames(ribos_avs) <- names(ribo_avs_ur)
ribo_proportions <- merge(mat, mat2, by = "row.names", all =TRUE)
row.names(ribo_proportions) <- ribo_proportions[,1]
ribo_proportions <- ribo_proportions[2:27]
write.csv(ribo_proportions, paste0(DATAPATH, "2017-02-17-proportions-of-ribosomal-gene-expression.csv", row.names = TRUE))

