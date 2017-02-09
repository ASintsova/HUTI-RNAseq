source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
biocLite("gage")
biocLite("org.EcK12.eg.db")
biocLite("AnnotationDbi")
library(AnnotationDbi)
library(pathview)
library(gage)
library(org.EcK12.eg.db)

library(readr)
DE_genes <- read.csv("~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/2017-01-23-core-genes-DE-edited.csv", header = TRUE)
de_genes <- DE_genes[,2:7]
row.names(de_genes)<- make.names(DE_genes[,8], TRUE)



#create a vector of gene_symbols
symbols <- row.names(de_genes)

head(symbols)

#run id2eg function 
ec.eg <- id2eg(ids = symbols, category = gene.idtype.list[1], org = "EcK12")

colnames(ec.eg) <- c("gene_name", "ENSMBL")
head(ec.eg)
n <- c(3:16)
gene_list <- cbind(ec.eg[,2], de_genes) # includes ensmbl number thing
colnames(gene_list) <- c("ENSMBL", colnames(de_genes))
head(gene_list)

tail(gene_list)


gene_list <- gene_list[complete.cases(gene_list),]

analyze_list <- as.matrix(gene_list[,3])
rownames(analyze_list) <- gene_list[,1]
colnames(analyze_list) <- c("LFC")

#analyze_list <- as.matrix(gene_list[,2:15])
#rownames(analyze_list)<- gene_list[,1]
#colnames(analyze_list)<- colnames(gene[,3:16])

#kegg annotation
kegg.ec <- kegg.gsets(species = "eco", id.type = "entrez")
kg.ec <- kegg.ec$kg.sets[kegg.ec$sigmet.idx]

#gage
path.analysis <- gage(analyze_list, gsets = kg.ec)

pathways_upregulated <- path.analysis$greater
pathways_downregulated <- path.analysis$less

write.csv(pathways_upregulated, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/2017-01-23-core-genes-pathways-UP.csv", row.names = TRUE)
write.csv(pathways_downregulated, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/2017-01-23-core-genes-pathways-DOWN.csv", row.names = TRUE)

#pathview
map.pathway <- as.matrix(gene[,n])
rownames(map.pathway) <- symbols
pv.out <- pathview(gene.data = map.pathway, pathway.id = "00230", species = "eco", gene.idtype = gene.idtype.list[1], same.layer = F)
