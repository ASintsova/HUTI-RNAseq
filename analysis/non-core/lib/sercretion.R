dat = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/data/"
fig = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/figures/"
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(pathview)
library(gage)
library(org.EcK12.eg.db)
library(readr)


gene_info <- read.table(paste0(dat, "non_core_gene_info.txt"), header = TRUE, row.names = 1, sep = "\t")
cmb.sig.genes <- read.csv(paste0(dat, "combined_non_core_de_genes.csv"), header = TRUE, row.names = 1)
non.core <- merge(cmb.sig.genes, gene_info[, c(2,4)], by = "row.names")
rownames(non.core)<- non.core[,1]
non.core <- non.core[,2:9]

core <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/DE-seq-results-all-core-edited.csv",
                 row.names = 1)


kegg.ec <- kegg.gsets(species = "eco", id.type = "eco")
kg.ec <- kegg.ec$kg.sets[kegg.ec$sigmet.idx]

#t2ss eco03070

t2ss <- kg.ec$`eco03070 Bacterial secretion system`

non.core.t2ss <- intersect(t2ss, non.core[,8])

core.t2ss <- intersect(t2ss, core[,9])


all.t2ss.core <- core[core$bnum %in% core.t2ss,]

all.t2ss.core <- all.t2ss.core[order(all.t2ss.core$log2FoldChange, decreasing = TRUE),]
write.csv(all.t2ss.core, paste0(dat, 'all-t2ss-core.csv'))
t2ss.core.names <- rownames(all.t2ss.core[all.t2ss.core$log2FoldChange > 0.78,])


t2ss.noncore.names <- rownames(non.core[non.core$bnum %in% non.core.t2ss,])

names <- c(t2ss.core.names, t2ss.noncore.names)

write.table(as.data.frame(names), paste0(dat, "secretion_gene_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

all.locs <- read.table(paste0(dat, "secretion_ed.gff"), sep = "\t")
get_id <- function(x){
  gi <- unlist(strsplit(unlist(strsplit(x, c(";")))[1], "="))[2]
  return (gi)
}

g.ids<- sapply(as.character(all.locs[,3]), get_id)
all.locs[,3]<- g.ids
colnames(all.locs)<- c("genome", "location", "gene")


hm56.t2ss <- all.locs[all.locs$genome == "HM56_1",]
hm56.t2ss.locs <- hm56.t2ss$location
names(hm56.t2ss.locs)<- hm56.t2ss$gene

hm56.t2ss.locs


hist(hm56.pc2.locs, 100, xlim = c(0, 4000000))
abline(v = 3177816, col = 'blue')
abline(v = 764444, col = 'red')


hm14.t2ss <- all.locs[all.locs$genome == "HM14_1",]
hm14.t2ss.locs <- hm14.t2ss$location
names(hm14.t2ss.locs)<- hm14.t2ss$gene

hm14.t2ss.locs


hist(hm14.t2ss.locs, 100, xlim = c(0, 4400000))
abline(v = 1465521, col = 'blue')
abline(v = 3831641, col = 'red')



