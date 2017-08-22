source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
biocLite("gage")
biocLite("org.EcK12.eg.db")
library(pathview)
library(gage)
library(org.EcK12.eg.db)


#read in file with log2FC
gene_symbol <- read.table("gene_symbols_fold_change.csv", header = TRUE, sep = ",")

#create a vector of gene_symbols
symbols <- gene_symbol[,1]
head(symbols)
#run id2eg function 
ec.eg <- id2eg(ids = symbols, category = gene.idtype.list[1], org = "EcK12")
gene_list <- matrix(c(as.numeric(ec.eg[,2]), gene_symbol[,2]), ncol = 2, byrow = FALSE)
head(gene_list)
gene_list <- gene_list[complete.cases(gene_list),]
dim(gene_list)
analyze_list <- as.matrix(gene_list[,2])
rownames(analyze_list) <- gene_list[,1]
colnames(analyze_list)<- "FC"
kegg.ec <- kegg.gsets(species = "eco", id.type = "entrez")
kg.ec <- kegg.ec$kg.sets[kegg.ec$sigmet.idx]
path.analysis <- gage(analyze_list, gsets = kg.ec)

#pv.out <- pathview(gene.data = sym2, pathway.id = "01501", species = "eco", gene.idtype = gene.idtype.list[1])
