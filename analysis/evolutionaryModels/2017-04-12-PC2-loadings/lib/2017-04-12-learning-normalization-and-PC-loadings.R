data ="/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-04-12-PC2-loadings/data/"

fig = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/2017-04-12-PC2-loadings/figures/"

#From Steve

#Inverse rank normalization

#From Steve
invnorm = function(x) { qnorm((rank(x,na.last="keep",ties.method="random")-0.5)/sum(!is.na(x))); }



#From internet

my.invnorm = function(x)
{
 res = rank(x)
 res = qnorm(res/(length(res)+0.5))
 return(res)
}

#if both denominator and numerator are large numbers this should not make a difference

#playing around with my data

#import normalized counts & expreiment_info

counts <- read.csv(paste0(data, "normalized_counts_core.csv"), row.names = 1)
info <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/data/Patient_meta_info.csv",
                row.names = 1)


 hist(counts[,1]) # very far from normal




 hist(log2(counts[,1] +1)) #much closer to normal

 test1 <- my.invnorm(counts[,1])
 hist(test1)#looks pretty damn normal

 test2 <-invnorm(counts[,1])
 hist(test2) #histogram looks best

 inv.norm.c <- matrix(NA, ncol = 26, nrow = 2502, byrow = FALSE)
 rownames(inv.norm.c)<- rownames(counts)
 colnames(inv.norm.c) <- colnames(counts)

 for (i in 1:ncol(counts)){
   x <- counts[,i]
   norm.x <- invnorm(x)
   inv.norm.c[,i] <- norm.x

 }

write.csv(inv.norm.c, paste0(data, "inverse_normalized_counts.csv"))


#Performing PCA analysis on these inverse normalzied counts
library(stats)
library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)



t.norm.counts <- t(inv.norm.c)
t.norm.counts <- merge(t.norm.counts, info[, c(2,4,5)], by ="row.names")
t.norm.counts <- t.norm.counts[,2:2506]
rownames(t.norm.counts)<- colnames(inv.norm.c)
media <- t.norm.counts[, 2503]
his <- t.norm.counts[,2504]
psg <- t.norm.counts[,2505]
all.pca <- prcomp(t.norm.counts[,1:2502])

#plotting all samples, labeling by media

g <- ggbiplot(all.pca, obs.scale = 1, var.scale = 1,
             groups = media, ellipse = TRUE, circle = FALSE,
             var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
                             legend.position = 'top')
ggsave(paste0(fig, "all_sample_pca_media.png"),g)

#plotting all samples by phylogroup

g <- ggbiplot(all.pca, obs.scale = 1, var.scale = 1,
              groups = psg, ellipse = TRUE, circle = FALSE,
                var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
                legend.position = 'top')
ggsave(paste0(fig, "all_sample_pca_psg.png"),g)


#urine
odd <- seq(1, 26, 2)
t.urine <- t.norm.counts[odd, ]
u.psg <- t.urine[, 2505]
u.his <- t.urine[, 2504]
u.pca <- prcomp(t.urine[,1:2502])

#plotting by psg
g <- ggbiplot(u.pca, obs.scale = 1, var.scale = 1,
              groups = u.psg, ellipse = TRUE, circle = FALSE,
                var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
                legend.position = 'top')
ggsave(paste0(fig, "urine_sample_pca_psg.png"),g)


#plotting by history
g <- ggbiplot(u.pca, obs.scale = 1, var.scale = 1,
              groups = u.his, ellipse = TRUE, circle = FALSE,
                var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
                legend.position = 'top')
ggsave(paste0(fig, "urine_sample_pca_his.png"),g)

#uti

even <- seq(2, 26, 2)
t.uti <- t.norm.counts[even, ]
uti.psg <- t.uti[, 2505]
uti.his <- t.uti[, 2504]
uti.pca <- prcomp(t.uti[,1:2502])

#plotting by psg
g <- ggbiplot(uti.pca, obs.scale = 1, var.scale = 1,
              groups = uti.psg, ellipse = TRUE, circle = FALSE,
                var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
                legend.position = 'top')
ggsave(paste0(fig, "uti_sample_pca_psg.png"),g)


#plotting by history
g <- ggbiplot(uti.pca, obs.scale = 1, var.scale = 1,
              groups = uti.his, ellipse = TRUE, circle = FALSE,
                var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
                legend.position = 'top')
ggsave(paste0(fig, "uti_sample_pca_his.png"),g)



#PCA loadings
annot = info[, c("MEDIA", "PRED_PHYLO")]
test2 <- prcomp(t(inv.norm.c))
test3 <- test2$rotation[,2][order(test2$rotation[,2], decreasing = TRUE)][1:10]#PC2
names <- names(test3)
load.genes <- inv.norm.c[names,]

write.csv(load.genes, paste0(data, "pc2_loading_genes.csv"))

pheatmap(load.genes, annotation_col = annot, cellheight = 10, cellwidth = 10)#,
           # filename = paste0(fig, 'pc2_loadings.png'))

test4<- test2$rotation[,1][order(test2$rotation[,1], decreasing = TRUE)][1:10]#PC1
names <- names(test4)
load.genes1 <- inv.norm.c[names,]

pheatmap(load.genes, annotation_col = annot, cellheight = 10, cellwidth = 10,
          filename = paste0(fig, "pc1_loadings.png" ))

write.csv(load.genes1, paste0(data, "pc1_loading_genes.csv"))


g.c <-log2(counts[names,]+1)
pheatmap(g.c, annotation_col = annot, cellheight = 10, cellwidth = 10)#,

