
library(phytools)
library(DESeq2)
library(pheatmap)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")


DATAPATH = "~/git_repos/HUTI-RNAseq/analysis/growth_rate/data/"

FIG = "~/git_repos/HUTI-RNAseq/analysis/growth_rate/figures/"


all_data <- LoadHTSeqCountsandTree()


PTR <- matrix(c(2.12, 2.29, 1.78, 2.38, 2.44, 1.79, 2.08, 1.21, 2.22, 2.32, 23.2, 20.0, 34.4, 18.7, 17.9), ncol = 3, nrow = 5, byrow = FALSE)
rownames(PTR) <- c("HM14", "HM17", "HM56", "HM66", "HM68")
colnames(PTR) <- c("PTR", "GR", "Doub.time")


counts <- all_data$uti_counts_l[, c("UTI14", "UTI17", "UTI56", "UTI66", "UTI68")]



row_sub = apply(counts, 1, function(row) all(row > 2 ))
counts <- counts[row_sub,]


ptr_genes <- matrix (NA, ncol=3, nrow=3000)
colnames(ptr_genes) <- c('gene', 'pval', 'int')
for (i in 1:nrow(counts)){
  m <- PTR
  m <- cbind(m, as.numeric(counts[i,  c("UTI14", "UTI17", "UTI56", "UTI66", "UTI68")]))
  colnames(m)[4] <- c("mean.exprs")
  m <- as.data.frame(m)
  mod <- lm(PTR ~ mean.exprs, data = m)
  z <- anova(mod)
  if (z$`Pr(>F)`[1] < 0.05){
    ptr_genes[i,1] <- rownames(counts[i,])
    ptr_genes[i,2]<- z$`Pr(>F)`[1]
    ptr_genes[i, 3] <- mod$coefficients[1]
    
  }
  
}
ptr_genes <- as.data.frame(ptr_genes)
ptr_genes <- ptr_genes[complete.cases(ptr_genes),]




#pos_genes <- ptr_genes[as.numeric(as.character(ptr_genes[,3])) > 0,]
rownames(ptr_genes) <- ptr_genes[,1]
ptr_genes <- ptr_genes[,c(2,3)]
ptr_genes <- ptr_genes[order(ptr_genes$pval),]
write.csv(ptr_genes, paste0(DATAPATH, "ptr_genes_correlated.csv"), row.names = TRUE)




ptr_gene_locs <- read.table(paste0(DATAPATH, "first_go_at_genome_locations.txt"), header = TRUE)

hm14.ptr.genes <- ptr_gene_locs[ptr_gene_locs$genome == "HM14",]
hm14.nhoA <- hm14.ptr.genes[hm14.ptr.genes$gene_id == "18964_nhoA", 4]
hm14.ptr.genes <- hm14.ptr.genes[!grepl("15607_acrB_2", hm14.ptr.genes$gene_id),]
rownames(hm14.ptr.genes) <- hm14.ptr.genes$gene_id
hm14.ptr.locs <- hm14.ptr.genes[,4]
names(hm14.ptr.locs) <-hm14.ptr.genes$gene_id

hist(hm14.ptr.locs, 500)
abline(v = 1465521, col = 'blue')
abline(v = 3831641, col = 'red')
abline(v = hm14.nhoA, col = 'green')

hm14.origin.genes <- hm14.ptr.locs[hm14.ptr.locs > 3831641 & hm14.ptr.locs < 4000000]

hm3.ptr.genes <- ptr_gene_locs[ptr_gene_locs$genome == "HM3",]

hm3.ptr.genes <- hm3.ptr.genes[!grepl("15607_acrB_2", hm3.ptr.genes$gene_id),]
rownames(hm3.ptr.genes) <- hm3.ptr.genes$gene_id
hm3.ptr.locs <- hm3.ptr.genes[,4]
names(hm3.ptr.locs) <-hm3.ptr.genes$gene_id

hist(hm3.ptr.locs, 500)
abline(v = 15050, col = 'blue')
abline(v = 2179454, col = 'red')
hm3.origin.genes <- hm14.ptr.locs[hm14.ptr.locs > 2179454 & hm14.ptr.locs < 2500000]



hm17.ptr.genes <- ptr_gene_locs[ptr_gene_locs$genome == "HM17",]

hm17.ptr.genes <- hm17.ptr.genes[!grepl("17383_intA_2", hm17.ptr.genes$gene_id),]
rownames(hm17.ptr.genes) <- hm17.ptr.genes$gene_id
hm17.ptr.locs <- hm17.ptr.genes[,4]
names(hm17.ptr.locs) <-hm17.ptr.genes$gene_id

hist(hm17.ptr.locs, 250)
abline(v = 5086096, col = 'blue')
abline(v = 2560407, col = 'red')
hm3.origin.genes <- hm14.ptr.locs[hm14.ptr.locs > 2179454 & hm14.ptr.locs < 2500000]








top_genes <- head(ptr_genes[order(ptr_genes$pval),])
names <- rownames(top_genes)

plotCor <- function(name){
  
  m <- cbind(PTR, as.numeric(counts[name,]))
  colnames(m)[4]<- c("gene.exprs")
  m <- as.data.frame(m)
  print (name)
  
  p <- ggplot(m, aes(x=PTR, y=gene.exprs))+
    geom_point(shape=1)+
    geom_smooth(method=lm)
  plot(p)
  return (p)
  
}


for (i in 1:length(names)){
  
  m <- cbind(PTR, as.numeric(counts[names[i],]))
  colnames(m)[4]<- c("gene.exprs")
  m <- as.data.frame(m)
  #print (names[i])
  
  p <- ggplot(m, aes(x=PTR, y=gene.exprs))+
    geom_point(shape=1)+
    geom_smooth(method=lm) + ggtitle(names[i])
  plot(p)
  #ggsave(paste0(FIG, names[i], "_vs_ptr.jpg"), p)
  Sys.sleep(4)
}

p <- ggplot(ptr_vs_ribo_uti, aes(x=PTR, y=mean.ribo.exprs))+
  geom_point(shape=1)+
  geom_smooth(method=lm)





