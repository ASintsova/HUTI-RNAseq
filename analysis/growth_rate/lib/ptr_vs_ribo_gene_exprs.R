#Calculating lambda on the same dataset

library(phytools)
library(DESeq2)
library(pheatmap)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")


DATAPATH = "~/git_repos/HUTI-RNAseq/analysis/DE/2017-02-17-ribosome-gene-expression/data/"
FIG = "~/git_repos/HUTI-RNAseq/analysis/growth_rate/figures/"

all_data <- LoadHTSeqCountsandTree()
ribo_genes <- as.vector(read.table(paste0(DATAPATH, "ribosomal_gene_ids.txt"))[,1])

counts <- all_data$counts
ribo_expr <- counts[ribo_genes,]

ribo_expr <- ribo_expr[complete.cases(ribo_expr),]
ribo_log <- log2(ribo_expr +1)
annot <- all_data$exp_info[, c("MEDIA", "PRED_PHYLO")]
ribo_log_mod <- ribo_log[!rownames(ribo_log) %in% c("15746_rplS", "15749_rpsP"),]#, "17881_ykgO"),]#, "18981_sra", "16752_rpmG",
# "17882_ykgM", "18114_ribosome-associated_protein" ),]
names <- row.names((ribo_log_mod))
disp_counts <- round(counts[names,]/1000, 1)

PTR <- matrix(c(2.11, 2.28, 1.77, 2.39, 2.43, 2.56, 2.97, 1.72, 3.23, 3.33, 16.3, 14.0, 24.1, 12.9, 12.5), ncol = 3, nrow = 5, byrow = FALSE)
rownames(PTR) <- c("HM14", "HM17", "HM56", "HM66", "HM68")
colnames(PTR) <- c("PTR", "GR", "Doub.time")

ptr_ribo_ur <- ribo_log_mod[, c("UR14", "UR17", "UR56", "UR66", "UR68")]

ptr_vs_ribo_ur <- cbind(PTR, colMeans(ptr_ribo_ur))

colnames(ptr_vs_ribo_ur)[4] <- c("mean.ribo.exprs")
ptr_vs_ribo_ur <- as.data.frame(ptr_vs_ribo_ur)

row_sub = apply(ribo_log_mod, 1, function(row) all(row > 5 ))

ribo.exprs <- ribo_log_mod[row_sub,]

ptr_ribo_ur <- ribo.exprs[, c("UR14", "UR17", "UR56", "UR66", "UR68")]

ptr_vs_ribo_ur <- cbind(PTR, colMeans(ptr_ribo_ur))
colnames(ptr_vs_ribo_ur)[4] <- c("mean.ribo.exprs")
ptr_vs_ribo_ur <- as.data.frame(ptr_vs_ribo_ur)
#plot(ptr_vs_ribo_ur$PTR, ptr_vs_ribo_ur$Doub.time)
#plot(ptr_vs_ribo_ur$Doub.time, ptr_vs_ribo_ur$mean.ribo.exprs)




ptr_ribo_uti <- ribo.exprs[, c("UTI14", "UTI17", "UTI56", "UTI66", "UTI68")]

ptr_vs_ribo_uti <- cbind(PTR, colMeans(ptr_ribo_uti))
colnames(ptr_vs_ribo_uti)[4] <- c("mean.ribo.exprs")
ptr_vs_ribo_uti <- as.data.frame(ptr_vs_ribo_uti)
#plot(ptr_vs_ribo_uti$PTR, ptr_vs_ribo_uti$Doub.time)
#plot(ptr_vs_ribo_uti$PTR, ptr_vs_ribo_uti$mean.ribo.exprs)


par(mfrow=c(1,2))

plot(ptr_vs_ribo_ur$Doub.time, ptr_vs_ribo_ur$mean.ribo.exprs)
plot(ptr_vs_ribo_uti$Doub.time, ptr_vs_ribo_uti$mean.ribo.exprs)


p <- ggplot(ptr_vs_ribo_uti, aes(x=PTR, y=mean.ribo.exprs))+
    geom_point(shape=1)+
    geom_smooth(method=lm)


mod <- lm (PTR ~ mean.ribo.exprs, ptr_vs_ribo_uti)
z <- anova(mod)

ggplot(ptr_vs_ribo_ur, aes(x=PTR, y=mean.ribo.exprs))+
  geom_point(shape=1)+
  geom_smooth(method=lm) + ggtitle('mean.urine.ribosomal.gene.expr')

test2 <- cbind(ptr_vs_ribo_uti, log10(ptr_vs_ribo_uti$Doub.time))
colnames(test2)[5] <- c("DT")
plot(test2$DT,log10(test2$mean.ribo.exprs))


for (i in 1:nrow(ribo_log_mod)){
  print (rownames(ribo_log_mod[i, ]))
  m <- cbind(PTR, as.numeric(ribo_log_mod[i, c("UTI14", "UTI17", "UTI56", "UTI66", "UTI68")]))
  colnames(m)[4] <- c("mean.ribo.exprs")
  m <- as.data.frame(m)
  print (m)
  plot(m$Doub.time, m$mean.ribo.exprs, main = rownames(ribo_log_mod[i,]))
  Sys.sleep(4)
  #mod <- lm(PTR ~ mean.ribo.exprs, data = m)
  
}




ptr_ribo <- ptr_vs_ribo_ur
ptr_ribo <-rbind(ptr_ribo, ptr_vs_ribo_uti)
ptr_ribo <- cbind(ptr_ribo, c(rep("UR", 5),rep("UTI",5)))
rownames(ptr_ribo)<- c("UR14", "UR17", "UR56", "UR66", "UR68","UTI14", "UTI17", "UTI56", "UTI66", "UTI68" )
colnames(ptr_ribo)[5] <- c("cond")
p <- ggplot(ptr_ribo, aes(x=PTR, y=mean.ribo.exprs, color=cond))+
  geom_point(shape=1)+
  scale_color_hue(l=50)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()

ggsave(paste0(FIG, "2017-03-28-mean-ribo-exprs-vs-PTR.jpg"),p )





counts <- all_data$uti_counts_l[, c("UTI14", "UTI17", "UTI56", "UTI66", "UTI68")]
row_sub = apply(counts, 1, function(row) all(row > 2 ))
counts <- counts[row_sub,]
ptr_genes <- matrix (NA, ncol=3, nrow=3000)
colnames(ptr_genes) <- c('gene', 'pval', 'int')
for (i in 1:nrow(counts)){
  m <- PTR
  m <- cbind(m, as.numeric(counts[i,  c("UTI14", "UTI17", "UTI56", "UTI66", "UTI68")]))
  colnames(m)[4] <- c("mean.ribo.exprs")
  m <- as.data.frame(m)
    mod <- lm(PTR ~ mean.ribo.exprs, data = m)
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
top_genes <- head(ptr_genes[order(ptr_genes$pval),])
names <- rownames(top_genes)

for (i in 1:length(names)){

  m <- cbind(PTR, as.numeric(counts[names[i],]))
  colnames(m)[4]<- c("gene.exprs")
  m <- as.data.frame(m)
  print (names[i])
  p <- ggplot(m, aes(x=PTR, y=gene.exprs))+
    geom_point(shape=1)+
    geom_smooth(method=lm)
  plot(p)
  ggsave(paste0(FIG, names[i], "_vs_ptr.jpg"), p)
  Sys.sleep(4)
}

p <- ggplot(ptr_vs_ribo_uti, aes(x=PTR, y=mean.ribo.exprs))+
  geom_point(shape=1)+
  geom_smooth(method=lm)




#2017-04-11


rib_locs <- read.table(paste0(DATAPATH, "ribosomal_go_at_genome_locations.txt"), header = TRUE)

hm17.rib.genes <- rib_locs[rib_locs$genome == "HM17",]

#hm17.rib.genes <- hm17.rib.genes[!grepl("17383_intA_2", hm17.ptr.genes$gene_id),]
rownames(hm17.rib.genes) <- hm17.rib.genes$gene_id
hm17.rib.locs <- hm17.rib.genes[,4]
names(hm17.rib.locs) <-hm17.rib.genes$gene_id

hist(hm17.rib.locs, 500)
abline(v = 5086096, col = 'blue')
abline(v = 2560407, col = 'red')
abline(v=4200000, col = 'grey')

hm14.rib.genes <- rib_locs[rib_locs$genome == "HM14",]

#hm17.rib.genes <- hm17.rib.genes[!grepl("17383_intA_2", hm17.ptr.genes$gene_id),]
rownames(hm14.rib.genes) <- hm14.rib.genes$gene_id
hm14.rib.locs <- hm14.rib.genes[,4]
names(hm14.rib.locs) <-hm14.rib.genes$gene_id

hist(hm14.rib.locs, 500)
abline(v = 1465521, col = 'blue')
abline(v = 3831641, col = 'red')



hm56.rib.genes <- rib_locs[rib_locs$genome == "HM56",]

#hm17.rib.genes <- hm17.rib.genes[!grepl("17383_intA_2", hm17.ptr.genes$gene_id),]
rownames(hm56.rib.genes) <- hm56.rib.genes$gene_id
hm56.rib.locs <- hm56.rib.genes[,4]
names(hm56.rib.locs) <-hm56.rib.genes$gene_id

hist(hm56.rib.locs, 500)
abline(v = 3177816, col = 'blue')
abline(v = 764444, col = 'red')


par(mfrow = c(1,3))





#2017-04-19
ptr_ribo_uti

ptr_genes <- matrix (NA, ncol=2, nrow=3000)
colnames(ptr_genes) <- c('gene', 'pval')
for (i in 1:nrow(ptr_ribo_uti)){
  m <- PTR
  m <- cbind(m, as.numeric(ptr_ribo_uti[i,  ]))
  colnames(m)[4] <- c("mean.ribo.exprs")
  m <- as.data.frame(m)
  mod <- lm(PTR ~ mean.ribo.exprs, data = m)
  z <- anova(mod)
  
  ptr_genes[i,1] <- rownames(ptr_ribo_uti[i,])
  ptr_genes[i,2]<- z$`Pr(>F)`[1]
    

}
ptr_genes <- as.data.frame(ptr_genes)
ptr_genes <- ptr_genes[complete.cases(ptr_genes),]
ptr_genes <- ptr_genes[order(ptr_genes$pval),]
