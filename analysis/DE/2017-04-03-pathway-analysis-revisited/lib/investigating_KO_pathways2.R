#set up
library(AnnotationDbi)
library(pathview)
library(gage)
library(org.EcK12.eg.db)
library(readr)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")
DATAPATH2 = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03/DE_genes_KO_ids.txt"
PATHINFO = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/path.txt"
OUTPATH = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/"


de_ko <- read.csv(DATAPATH2, row.names = 1, header = TRUE)
cc_de_ko <- de_ko[complete.cases(de_ko),]

ko <- cc_de_ko[!duplicated(cc_de_ko[,1]),]

all_data <- LoadHTSeqCountsandTree()

#DE results from 2017-01-30
DE_genes = read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-01-30-htseq-counts/data/core_DEseq_results.csv",
                    header = TRUE, row.names = 1)
names <- rownames(DE_genes)
gene_names <-c()
for (name in names){
  gene_names <- append(gene_names, unlist(strsplit(name , "_"))[2])
}

de_genes<- cbind(DE_genes, gene_names)

data <- merge(ko, de_genes, by = "row.names")
data2 <- data[,c(2, 3,4,5,6,7,8,15)]
rownames(data2)<- data[,1]
colnames(data2) <- c("KO", "baseMean", "LFC", "lfcSE", "stat", "pval", "padj", "gene_names")
write.csv(data2, paste0(OUTPATH, "DE_genes_withKO_and_gene_names.csv"))


t1 <- data2[,c(1,3)]
t1$LFC <- as.numeric(as.character(t1$LFC))
t1$LFC <- as.numeric(as.character(t1$LFC)
rownames(t1)<- t1[,1]
t1 <- as.matrix(t1)

t1$LFC <- as.numeric(as.character(t1$LFC))
pv.out <- pathview(gene.data = map.pathway, pathway.id = "00230",  species = "eco",
                   gene.idtype = "kegg", outsuffix = "test", same.layer=F)







#Alanine, Aspartate, Glutamate 00250
path <- read.table(PATHINFO)
genes <- intersect(path[,1], data2[,1])

de_genes <- subset(data2, KO %in% genes)
names <- de_genes[,c(1, 8)]
names <- names[!rownames(names) %in% c("18677_purB"),]

DE_expr <- log2(all_data$counts[rownames(names),]+1)  #norm counts for all the diff expr genes

de_expr <- as.matrix(DE_expr)
disp_counts <- round(all_data$counts[rownames(names),]/100, 1)
rownames(de_expr) <- names[,2]



#Ala/Asp/Glu

eco00250 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00250, paste0(OUTPATH, "eco00250"))
#eco00250 <- read.csv(paste0(OUTPATH, "eco00250"), row.names = 1)
eco00250_mean <- matrix(colMeans(eco00250), ncol =26, nrow =1)

annot = all_data$exp_info[, c("MEDIA", "PRED_PHYLO")]
pheatmap(eco00250, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
          cutree_cols = 2 )




#Gly/Ser/Thr

eco00260 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00260, paste0(OUTPATH, "eco00260"))
eco00260_mean <- matrix(colMeans(eco00260), ncol =26, nrow =1)
pheatmap(eco00260, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )


#Val/Leu/Ileu

eco00290 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00290, paste0(OUTPATH, "eco00290"))
eco00290_mean <- matrix(colMeans(eco00290), ncol =26, nrow =1)


pheatmap(eco00290, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )



#Arg

eco00220 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00220, paste0(OUTPATH, "eco00220"))
eco00220_mean <- matrix(colMeans(eco00220), ncol =26, nrow =1)

pheatmap(eco00220, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )


#Cys/Met

eco00270 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00270, paste0(OUTPATH, "eco00270"))
eco00270_mean <- matrix(colMeans(eco00270), ncol =26, nrow =1)

pheatmap(eco00270, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )





aa <- rbind( eco00250_mean, eco00260_mean,eco00220_mean, eco00270_mean, eco00290_mean)

colnames(aa) <- colnames(eco00220)
rownames(aa) <- c("Ala/Asp/Glu", "Gly/Ser/Thr", "Arg", "Cys/Met", "Val/Leu/Ileu")
#write.csv(aa, paste0(OUTPATH, "all_aa"))
pheatmap(aa, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE, cluster_rows = FALSE,
         filename = paste0(OUTPATH, 'aa_acid_metabolism.jpg'))


#tca

eco00020 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00020, paste0(OUTPATH, "eco00020"))
eco00020_mean <- matrix(colMeans(eco00020), ncol =26, nrow =1)
pheatmap(eco00020, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )


#galactose
eco00052 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00052, paste0(OUTPATH, "eco00052"))
eco00052_mean <- matrix(colMeans(eco00052), ncol =26, nrow =1)

pheatmap(eco00052, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )

#fructose and mannose

eco00051 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00051, paste0(OUTPATH, "eco00051"))
eco00051_mean <- matrix(colMeans(eco00051), ncol =26, nrow =1)

pheatmap(eco00051, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )

#pentose phosphate pathway

eco00030 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00030, paste0(OUTPATH, "eco00030"))
eco00030_mean <- matrix(colMeans(eco00030), ncol =26, nrow =1)

pheatmap(eco00030, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )

carb_met <- rbind( eco00020_mean, eco00052_mean,eco00051_mean, eco00030_mean)
colnames(carb_met) <- colnames(eco00020)
rownames(carb_met) <- c("TCA", "Gal", "Fru/Man", "PentPhosph")


#purine



eco00230 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00230, paste0(OUTPATH, "eco00230"))
eco00230_mean <- matrix(colMeans(eco00230), ncol =26, nrow =1)
pheatmap(eco00230, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )

#pyrimidine

eco00240 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco00240, paste0(OUTPATH, "eco00240"))
eco00240_mean <- matrix(colMeans(eco00240), ncol =26, nrow =1)
pheatmap(eco00240, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2, filename = paste0(OUTPATH, "pyrimidine_uti.jpg") )



#NucMet

nuc_met <- rbind(eco00230_mean, eco00240_mean)
colnames(nuc_met) <- colnames(eco00230)
rownames(nuc_met) <- c("Purine", "Pyr")

#Translation

#ribo

eco03010 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco03010, paste0(OUTPATH, "eco03010"))
eco03010_mean <- matrix(colMeans(eco03010), ncol =26, nrow =1)
pheatmap(eco03010, annotation_col = annot, 
         cellwidth = 7, cellheight = 7,  cluster_cols =FALSE,
         cutree_cols = 2 )



transl <- eco03010_mean
colnames(transl) <- colnames(eco03010)
rownames(transl)<- c("Ribosome")


#DNA replication



eco03030 <- de_expr[, order(colnames(de_expr))]
write.csv(eco03030, paste0(OUTPATH, "eco03030"))
eco03030_mean <- matrix(colMeans(eco03030), ncol =26, nrow =1)
pheatmap(eco03030, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )
#NADA

##RNA poly


eco03020 <- de_expr[, order(colnames(de_expr))]
write.csv(eco03020, paste0(OUTPATH, "eco03020"))
ecoeco03020_mean <- matrix(colMeans(eco03020), ncol =26, nrow =1)

pheatmap(eco03020, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )



#tRNA

eco00970 <- de_expr[, order(colnames(de_expr))]
write.csv(eco00970, paste0(OUTPATH, "eco00970"))
eco00970_mean <- matrix(colMeans(eco00970), ncol =26, nrow =1)
pheatmap(eco00970, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )




#RNA degradation



eco03018 <- de_expr[, order(colnames(de_expr))]
write.csv(eco03018, paste0(OUTPATH, "eco03018"))
eco03018_mean <- matrix(colMeans(eco03018), ncol =26, nrow =1)

pheatmap(eco03018, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )



#Quorum sensing



eco02024 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco02024, paste0(OUTPATH, "eco02024"))
eco02024_mean <- matrix(colMeans(eco02024), ncol =26, nrow =1)
pheatmap(eco02024, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )

#Biofilm formation

eco02026 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco02026, paste0(OUTPATH, "eco02026"))
#eco02026 <- read.csv(paste0(OUTPATH, "eco02026"), row.names = 1)
eco02026_mean <- matrix(colMeans(eco00250), ncol =26, nrow =1)
pheatmap(eco02026, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2, filename = paste0(OUTPATH, 'biofilm_uti.jpg') )

#Cell community

cell_comm <- rbind(eco02024_mean, eco02026_mean)
colnames(cell_comm)<- colnames(eco02024)
rownames(cell_comm) <- c("Quorum", "Biofilm")



####METABOLISM

MET <- rbind(aa, carb_met, nuc_met, transl, cell_comm)

pheatmap(MET, annotation_col = annot, gaps_row = c(5, 9,11,12),
         gaps_col = 13,
         cellwidth = 15, cellheight = 15,  cluster_cols =FALSE,
         cluster_rows = FALSE,
         cutree_cols = 2, annotation_legend = FALSE,
         filename = paste0(OUTPATH, "metabolism_during_uti.jpg"))




#Bacterial Chemotaxis

eco02030 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco02030, paste0(OUTPATH, "eco02030"))
eco02030_mean <- matrix(colMeans(eco02030), ncol =26, nrow =1)
pheatmap(eco02030, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE,
         cutree_cols = 2 )




#Transporters


eco02010 <- de_expr[, order(colnames(de_expr))]
#write.csv(eco02010, paste0(OUTPATH, "eco02010"))
eco02010_mean <- matrix(colMeans(eco02010), ncol =26, nrow =1)
pheatmap(eco02010, annotation_col = annot, 
         cellwidth = 10, cellheight = 10,  cluster_cols =FALSE, annotation_legend = FALSE,
         cutree_cols = 2, filename = paste0(OUTPATH, "abc_transp_uti.jpg"))



