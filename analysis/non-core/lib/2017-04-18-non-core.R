
#Set up

dat = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/data/"
fig = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/figures/"
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(pathview)
library(gage)
library(org.EcK12.eg.db)
library(readr)

#Load in DESeq counts
counts <- read.csv(paste0(dat, "2017-02-01-combined_counts.csv"), row.names = 1)
strain_names <- c("UR1", "UTI1", "UR2", "UTI2", "UR3", "UTI3", "UR4", "UTI4", "UR5", "UTI5", "UR6", "UTI6", "UR7", "UTI7",
  "UR8", "UTI8", "UR9", "UTI9", "UR10", "UTI10", "UR11", "UTI11", "UR12", "UTI12")

exp_info <- read.csv(paste0(dat, "non-core-exp-info.csv"), row.names = 1)
cmb.sig.genes <- data.frame()


####

#Trying to make this into an actual script

genome.num <- 5:12


for (i in genome.num){
  
  print(i)
  num.missing <- 26 - i*2
  print(num.missing)
  all.count <- counts[rowSums(is.na(counts)) == num.missing,]
  cts <- matrix(NA, ncol = i*2, nrow = nrow(all.count))
  rownames(cts)<-rownames(all.count)
  colnames(cts)<- strain_names[1:ncol(cts)]
  for (j in 1:nrow(all.count)){
    mod <- all.count[j,]
    
    cts[j, ]<- mod[!is.na(mod)]
    
  }
  write.csv(cts, paste0(dat, "counts_", i, "_genomes.csv"))
  
  info <- exp_info[rownames(exp_info) %in% colnames(cts),,drop = FALSE]
  
  cts.dds <- DESeqDataSetFromMatrix(countData = cts, colData = info,
                                    design = ~MEDIA)
  
  cts.dds.transformed <- rlog(cts.dds, blind=FALSE)
  pcaData <- plotPCA(cts.dds.transformed, intgroup = c("MEDIA"), returnData = TRUE)
  percentVar <- round(100*attr(pcaData, "percentVar"))
  g <- ggplot(data = pcaData, aes(x=PC1, y=PC2, color = MEDIA)) + 
    geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% var")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  ggsave(paste0(fig, "pca_", i, "_genomes.png"), g)
  
  cts.de <- DESeq(cts.dds)
  
  cts.sign <- subset(results(cts.de, lfcThreshold = 1), padj < 0.05)
  cts.sign <- cts.sign[order(cts.sign$log2FoldChange, decreasing = TRUE),]
  write.csv(cts.sign, paste0(dat, "sig_de_", i, "_genomes.csv"))
  
  cmb.sig.genes <- rbind(cmb.sig.genes, cts.sign)
  
  
}


cmb.sig.genes <- cmb.sig.genes[order(cmb.sig.genes$log2FoldChange, decreasing = TRUE),]
write.csv(cmb.sig.genes, paste0(dat, "combined_non_core_de_genes.csv"))

gene_info <- read.table(paste0(dat, "non_core_gene_info.txt"), header = TRUE, row.names = 1, sep = "\t")


non.core <- merge(cmb.sig.genes, gene_info[, c(2,4)], by = "row.names")
rownames(non.core)<- non.core[,1]
non.core <- non.core[,2:9]

kegg.ec <- kegg.gsets(species = "eco", id.type = "eco")
kg.ec <- kegg.ec$kg.sets[kegg.ec$sigmet.idx]

p.a <- non.core[,2]
names(p.a) <- non.core[,8]

path.analysis <- gage(p.a, gsets = kg.ec)

head(path.analysis$greater[,1:5], 10)
head(path.analysis$less[,1:5], 10)


pv.out <- pathview(gene.data = p.a, pathway.id = "03070", species = "eco", gene.idtype = "ENTREZID", same.layer = F)

pv.out <- pathview(gene.data = p.a, pathway.id = "02010", species = "eco", gene.idtype = "ENTREZID", same.layer = F)



#want to combine core and non-core

core <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/DE_genes_edited.csv",
                 row.names = 1)

all.genes <- rbind(core[,c(2,7,9)], non.core[,c(2,7,8)])

write.csv(all.genes, paste0(dat, "all_de_genes.csv"))

p.a <- all.genes[,1]
names(p.a) <- all.genes[,2]

path.analysis <- gage(p.a, gsets = kg.ec)

head(path.analysis$greater[,1:5], 10)
head(path.analysis$less[,1:5], 10)

pv.out <- pathview(gene.data = p.a, pathway.id = "00052", species = "eco", gene.idtype = "ENTREZID", same.layer = F)


core.p.a <- core[,2]
names(core.p.a)<-core[,9]
pv.out <- pathview(gene.data = core.p.a, pathway.id = "02010", species = "eco", gene.idtype = "ENTREZID", same.layer = F)






