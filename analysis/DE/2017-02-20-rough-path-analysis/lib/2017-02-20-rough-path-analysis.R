# Since I still don't have identifiers for each gene going to have to continue working with 
#gene names
DATAPATH ="/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-02-20-rough-path-analysis/data/"
FIGUREPATH = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-02-20-rough-path-analysis/figures/"


#set up
library(AnnotationDbi)
library(pathview)
library(gage)
library(org.EcK12.eg.db)
library(readr)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")

all_data <- LoadHTSeqCountsandTree()

#DE results from 2017-01-30
DE_genes = read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-01-30-htseq-counts/data/core_DEseq_results.csv",
                    header = TRUE, row.names = 1)

#added rough gene names
names <- rownames(DE_genes)
gene_names <-c()
for (name in names){
        gene_names <- append(gene_names, unlist(strsplit(name , "_"))[2])
}

de_genes<- cbind(DE_genes, gene_names)
write.csv(de_genes, paste0(DATAPATH,"2017-02-20-DE-genes-with-gene-names.csv"), row.names = TRUE)

#pathway analysis

#run id2eg function 
ec.eg <-id2eg(ids = gene_names, category = gene.idtype.list[1], org = "EcK12")
#Note: 82 of 559 unique input IDs unmapped
head(ec.eg)
#n <- c(3:16)
de_genes <- cbind(ec.eg[,2], de_genes) # includes ensmbl number thing
colnames(de_genes) <- c("ENTREZID", colnames(de_genes)[2:8])
head(de_genes)

tail(de_genes)
gene_list <- de_genes[complete.cases(de_genes),]

DE_names <- row.names(gene_list)
DE_expr <- log2(all_data$counts[DE_names,]+1)  #norm counts for all the diff expr genes
de_expr <- as.matrix(DE_expr)
rownames(de_expr) <- gene_list[,1]
cn = colnames(de_expr)
ur = grep ("UR", cn, ignore.case = TRUE)
uti = grep("UTI", cn, ignore.case = TRUE)

#kegg annotation
kegg.ec <- kegg.gsets(species = "eco", id.type = "entrez")
kg.ec <- kegg.ec$kg.sets[kegg.ec$sigmet.idx]

#gage
path.analysis <- gage(de_expr, gsets = kg.ec, ref = ur, samp = uti)
str(path.analysis, strict.width = 'wrap')
head(path.analysis$greater[,1:5], 4)
head(path.analysis$less[,1:5], 4)

pathways_upregulated <- path.analysis$greater
pathways_downregulated <- path.analysis$less

write.csv(pathways_upregulated, paste0(DATAPATH, "2017-02-20-pathways-upregulated.csv"), row.names = TRUE)
write.csv(pathways_downregulated, paste0(DATAPATH, "2017-02-20-pathways-downregulated.csv"), row.names = TRUE)

#Run the same with same.dir=F to look for groups perturbed both ways

path.analysis.no.direction <- gage(de_expr, gsets = kg.ec, ref = ur, 
                                   samp = uti, same.dir = F)
write.csv(path.analysis.no.direction$greater, 
          paste0(DATAPATH, "2017-02-20-pathways-distrubed-both-ways.csv"))

head(path.analysis.no.direction$greater[,1:5], 4)

#for small gene sets non-parametric test might be better:
path.analysis.MannWhit <- gage(de_expr, gsets = kg.ec, ref = ur, 
                               samp = uti, rank.test = T)
head(path.analysis.no.direction$greater[,1:5], 4)
#did not appear to change much in this case

#significantly changed gene_sets
path.analysis.sign <- sigGeneSet(path.analysis, outname=paste0(DATAPATH,"path.analysis.sign"))
write.csv(path.analysis.sign$greater, paste0(DATAPATH, "2017-02-20-pathways-upregulated-sign.csv"), row.names = TRUE)
write.csv(path.analysis.sign$less, paste0(DATAPATH, "2017-02-20-pathways-downregulated.csv"), row.names = TRUE)


#significant in 2 directions

path.analysis.no.direction.sign <- sigGeneSet(path.analysis.no.direction, outname=paste0(DATAPATH,"path.analysis.no.direction.sign"))
write.csv(path.analysis.no.direction.sign$greater, paste0(DATAPATH, "2017-02-20-pathways-two-directions-sign.csv"), row.names = TRUE)
path.analysis.no.direction.sign$greater

#essential gene sets

p.a.geneSets.up <- esset.grp(path.analysis$greater, de_expr, gsets = kg.ec, ref = ur,
                             samp = uti, test4up = T, output = T,
                             outname = paste0(DATAPATH, "2017-02-20-essentialGeneSets"),
                             make.plot = F)
                            

p.a.geneSets.down <- esset.grp(path.analysis$less, de_expr, gsets = kg.ec, ref = ur,
                             samp = uti, test4up = T, output = T,
                             outname = paste0(DATAPATH, "2017-02-20-essentialGeneSets-down"),
                             make.plot = F)

#something is not working out properly for this :S



#Pyrimidine metabolism: eco00240 
#Purine metabolism: eco00230
#Biosynthesis of amino acids: eco01230


#Butanoate metabolism: eco00650
#Galactose metabolism: eco00052
#Cysteine metabolism: eco00270
#Pyruvate: eco00620


#pathview
map.pathway <- as.matrix(gene[,n])
rownames(map.pathway) <- symbols
pv.out <- pathview(gene.data = map.pathway, pathway.id = "00230", species = "eco", gene.idtype = gene.idtype.list[1], same.layer = F)


