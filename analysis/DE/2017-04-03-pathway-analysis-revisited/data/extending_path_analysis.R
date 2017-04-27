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

de_counts <- all_data$counts[DE_names,]
de_counts <- as.matrix(de_counts)
rownames(de_counts) <- gene_list[,1]

cn = colnames(de_expr)
ur = grep ("UR", cn, ignore.case = TRUE)
uti = grep("UTI", cn, ignore.case = TRUE)

#kegg annotation
kegg.ec <- kegg.gsets(species = "eco", id.type = "entrez")

#signalling and metabolism

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

##For pathway visulization

d <- de_genes[,c(1, 2,3,7)]
d <- d[!duplicated(d[,1]),]
d <- d[complete.cases(d),]

rownames(d)<- as.character(d[,1])
#Purine metabolism: eco00230
purine <- p.a.geneSets.up$coreGeneSets$`eco00230 Purine metabolism`
pur_names <- eg2id(eg = purine, category = gene.idtype.list[1], org = "EcK12")
pur_expr <- de_expr[purine,]
test <- as.matrix(d[,3])
rownames(test) <- d[,1]
pv.out <- pathview(gene.data = map.pathway, pathway.id = "00230",  species = "eco",
                   gene.idtype = "kegg", outsuffix = "test", same.layer=F)


rownames(pur_expr)<- as.data.frame(pur_names)$SYMBOL
annot = all_data$exp_info[, c("MEDIA", "PRED_PHYLO")]
pheatmap(pur_expr, annotation_col = annot, cellheight = 20, cellwidth = 20)






#Butanoate metabolism: eco00650

#Cysteine metabolism: eco00270



#pathview
map.pathway <- as.matrix(gene[,n])
rownames(map.pathway) <- symbols
pv.out <- pathview(gene.data = map.pathway, pathway.id = "000270", species = "eco", gene.idtype = gene.idtype.list[1], same.layer = F)

#function to get genes from a pathway from your dataset, assumes existence of kg.ec, and de_expr

getPathwayAssociatedGenes <- function(path.assoc.genes){
        
        path.assoc.gene.names <- intersect(path.assoc.genes, rownames(de_expr))
        gene.expr <- de_expr[path.assoc.gene.names,]
        disp_counts <- round(de_counts[path.assoc.gene.names,]/1000, 1)
        rownames(gene.expr)<- as.data.frame(eg2id(eg = path.assoc.gene.names, category = gene.idtype.list[1], org = "EcK12"))$SYMBOL
        out <- list("gene_expression" = gene.expr, "counts" = disp_counts)
        return (out)
}



###Looking for genes in my dataself belonging to certain pathway
#Purine metabolism eco00230

pur_met <- kg.ec$eco00230
pur_expr_names <- intersect(pur_met, rownames(de_expr))
pur_expr <- de_expr[pur_expr_names,]
disp_counts <- round(de_counts[pur_expr_names,]/1000, 1)
rownames(pur_expr)<- as.data.frame(eg2id(eg = pur_expr_names, category = gene.idtype.list[1], org = "EcK12"))$SYMBOL
pheatmap(pur_expr, annotation_col = annot, display_numbers = disp_counts,
         cellwidth = 20, cellheight = 20,
         cutree_rows = 3, cutree_cols = 4, filename = paste0(FIGUREPATH, "2017-02-21-purine-metabolism.png"))







#Biosynthesis of amino acids: eco01230
aa.bio <- kg.ec$eco01230
aa.bio.data <- getPathwayAssociatedGenes(aa.bio)
aa.bio.expr <- aa.bio.data$gene_expression
aa.bio.counts <- aa.bio.data$counts
pheatmap(aa.bio.expr, annotation_col = annot, display_numbers = aa.bio.counts,
          cutree_cols = 2)

#Biosynthesis of cystein : kg.ec$`eco00270 Cysteine and methionine metabolism`
x <- kg.ec$`eco00270 Cysteine and methionine metabolism`
x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts)


#Pyrimidine metabolism: eco00240 
x <- kg.ec$`eco00240 Pyrimidine metabolism`
x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20, filename = paste0(FIGUREPATH, "2017-02-21-pyrimidine-metabolism.png"))

#Pyruvate: eco00620
x <- kg.ec$`eco00620 Pyruvate metabolism`
x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20)
         #filename = paste0(FIGUREPATH, "2017-02-21-pyrimidine-metabolism.png"))

#Galactose metabolism: eco00052

x <- kg.ec$`eco00052 Galactose metabolism`
x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20)


#Butanoate metabolism: eco00650

x <- kg.ec$`eco00650 Butanoate metabolism`
x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20)


#DNA replication

x <- kg.ec$`eco03030 DNA replication`
x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20)

# doesn't look like these have anything to do with DNA replication :S


#Fructose and mannose metabolism

x <-kg.ec$`eco00051 Fructose and mannose metabolism`

x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20)
# fructose genes are up

#Mismatch repair 
x <-kg.ec$`eco03430 Mismatch repair`

x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20)
#no sign differences dnaX is up, exoX is down

#Homologous recombination:
x <-kg.ec$`eco03440 Homologous recombination`

x.data <- getPathwayAssociatedGenes(x)
x.expr <- x.data$gene_expression
x.counts <- x.data$counts
pheatmap(x.expr, annotation_col = annot, display_numbers = x.counts,
         cellheight = 20, cellwidth = 20)

#not showing much, but it looks like I might be missing some genes, i.e. cannot finde
#recA in my dataset at all. 


#should look at different aa seperately: down: isoleucine, cys, arg, glt
         #filename = paste0(FIGUREPATH, "2017-02-21-purine-metabolism.png"))

#signaling 

kg.ec.sig <- kegg.ec$kg.sets[kegg.ec$sig.idx]

#gage
path.analysis.sig <- gage(de_expr, gsets = kg.ec.sig, ref = ur, samp = uti)
str(path.analysis.sig, strict.width = 'wrap')
head(path.analysis.sig$greater[,1:5], 10)
head(path.analysis.sig$less[,1:5], 10)


