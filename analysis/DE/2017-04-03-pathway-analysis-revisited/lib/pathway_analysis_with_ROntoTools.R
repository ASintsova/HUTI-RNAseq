source("https://bioconductor.org/biocLite.R")
biocLite("ROntoTools")

library(ROntoTools)

require(graph)
require(ROntoTools)
kpg <- keggPathwayGraphs("eco", relPercThresh = 0.0, updateCache = TRUE, verbose = TRUE)


kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
                      edgeWeightByType = list(activation = 1, inhibition = -1,
                                                expression = 1, repression = -1),
                      defaultWeight = 0)

kpn <- keggPathwayNames("eco")
core <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/DE-seq-results-all-core-edited.csv", row.names = 1, header = TRUE)
fc <- core$log2FoldChange[core$padj <= 0.01]
names(fc)<- paste0("eco:", core$bnum[core$padj <= 0.01])
pv <-core$padj[core$padj <= 0.01]
names(pv) <- paste0("eco:", core$bnum[core$padj <= 0.01])



kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
head(nodeWeights(kpg[["path:eco01501"]]))
ref <- paste0("eco:",core$bnum)
peRes <- pe(x = fc, graphs = kpg, ref = ref,  nboot = 200, verbose = TRUE)


p <- peRes@pathways[["path:eco02024"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)



peRes <- pe(x = fc2, graphs = kpg, ref = ref,  nboot = 1000, verbose = FALSE)
Summary(peRes, pathNames = kpn,  totalAcc = FALSE, totalPert = FALSE,
        pAcc = FALSE, order.by = "pPert")





#Going back to DESeq

library(DESeq2)
library (pheatmap)
library(RColorBrewer)
library (ggplot2)

#Importing the counts
combined_counts <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/2017-01-30-htseq-counts/data/2017-02-01-combined_counts.csv", header = TRUE, row.names = 1)

#Importing metadata
experiment_info <- read.csv("~/git_repos/HUTI-RNAseq/analysis/DE/data/Patient_meta_info.csv", header = TRUE, row.names =1)

head(combined_counts)

#Isolating Core Genes

core_genes = combined_counts[complete.cases(combined_counts),]
summary(core_genes)
(even <- seq(2, 26, 2))
(odd <- seq(1, 26, 2))


# Construct DESeqDataSet object from the matrix of counts and the experiment_info
(core_dds <- DESeqDataSetFromMatrix(countData = as.matrix(core_genes), colData = experiment_info,
                                    design = ~MEDIA))





#Data transform: DESeq2 specific, there are disadvantages to doing clustering analysis on RPKMs,
# and even normal log transforms
core_dds_transformed <- rlog(core_dds, blind=FALSE)#FALSE settin gwill means differences due to 
# design are not taken into account





# Look at other clustering methods: hclust??

plotPCA(core_dds_transformed, intgroup = c("MEDIA"))



core_DE <- DESeq(core_dds)


res <- results(core_DE) # this gives me lfc for every gene, going to use this for pathway analysis
#except that i first have to run it through python code

write.csv(res, "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/DE-seq-results-all-core.csv", row.names = TRUE)

core <- read.csv("../data/DE-seq-results-all-core-edited.csv", row.names = 1, header = TRUE)
fc2 <- core$log2FoldChange[core$padj <= 0.01]
names(fc2)<- paste0("eco:", core$bnum[core$padj <= 0.01])
pv <-core$padj[core$padj <= 0.01]
names(pv) <- paste0("eco:", core$bnum[core$padj <= 0.01])

kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
ref <- paste0("eco:",core$bnum)

peRes <- pe(x = fc2, graphs = kpg, ref = ref,  nboot = 500, verbose = FALSE)



results <- Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
        pAcc = FALSE,  comb.pv = NULL, order.by = "pPert")
write.csv(results,"/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/ROntoTools_output.csv")
kpg <- keggPathwayGraphs("eco", updateCache = TRUE, verbose = TRUE)



#kpg loading thing is not working for some reason


##Going back to gage



#set up
library(AnnotationDbi)
library(pathview)
library(gage)
library(org.EcK12.eg.db)
library(readr)
library(phytools)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")

all_data <- LoadHTSeqCountsandTree()

#DE results from 2017-01-30
DE_genes = read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/DE_genes_edited.csv",
                    header = TRUE, row.names = 1)

#added rough gene names
names <- rownames(DE_genes)

counts <- cbind(all_data$counts[names,], DE_genes$bnum)
c <- counts[complete.cases(counts),]
c<- c[!duplicated(c$`DE_genes$bnum`),]
rownames(c)<- c$`DE_genes$bnum`
c <- c[,c(1:26)]
#pathway analysis
de_expr <- log2(c + 1)
de_expr <- as.matrix(de_expr)

#de_counts <- all_data$counts[DE_names,]
#de_counts <- as.matrix(de_counts)
#rownames(de_counts) <- gene_list[,1]

cn = colnames(de_expr)
ur = grep ("UR", cn, ignore.case = TRUE)
uti = grep("UTI", cn, ignore.case = TRUE)

#kegg annotation
kegg.ec <- kegg.gsets(species = "eco")
kg.ec <- kegg.ec$kg.sets[kegg.ec$sigmet.idx]
kg.dise <- kegg.ec$kg.sets[kegg.ec$dise.idx]
kg.si <- kegg.ec$kg.sets[kegg.ec$sig.idx]
#gage
path.analysis <- gage(de_expr, gsets = kg.ec, ref = ur, samp = uti)

path.analysis.si <- gage(de_expr, gsets = kg.si, ref = ur, samp = uti)
str(path.analysis, strict.width = 'wrap')
head(path.analysis.si$greater[,1:5], 10)
head(path.analysis.si$less[,1:5], 10)

pathways_upregulated <- path.analysis$greater
pathways_downregulated <- path.analysis$less

write.csv(pathways_upregulated, "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/2017-04-06-pathways-upregulated.csv", row.names = TRUE)
write.csv(pathways_downregulated, "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/2017-04-06-pathways-downregulated.csv", row.names = TRUE)

#Run the same with same.dir=F to look for groups perturbed both ways

path.analysis.no.direction <- gage(de_expr, gsets = kg.ec, ref = ur, 
                                   samp = uti, same.dir = F)
write.csv(path.analysis.no.direction$greater, 
          "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/2017-04-06-pathways-distrubed-both-ways.csv")

head(path.analysis.no.direction$greater[,1:5], 10)


#significantly changed gene_sets
path.analysis.sign <- sigGeneSet(path.analysis, outname="/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/path.analysis.sign")
write.csv(path.analysis.sign$greater, "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/2017-04-06-pathways-upregulated-sign.csv", row.names = TRUE)
write.csv(path.analysis.sign$less, "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/2017-04-06-pathways-downregulated.csv", row.names = TRUE)


#essential gene sets

p.a.geneSets.up <- esset.grp(path.analysis$greater, de_expr, gsets = kg.ec, ref = ur,
                             samp = uti, test4up = T, output = T,
                             #outname = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2017-04-03-pathway-analysis-revisited/data/2017-04-06-essentialGeneSets",
                             make.plot = T)


p.a.geneSets.down <- esset.grp(path.analysis$less, de_expr, gsets = kg.ec, ref = ur,
                               samp = uti, test4up = T, output = T,
                               #outname = paste0(DATAPATH, "2017-02-20-essentialGeneSets-down"),
                               make.plot = F)

#pathview
map.pathway <- DE_genes$log2FoldChange
names(map.pathway)<- DE_genes$bnum

pv.out <- pathview(gene.data = map.pathway, pathway.id = "00240", gene.idtype = "ENTREZID", species = "eco", same.layer = F)

names <- names(uti_means[uti_means[,1] > 10,])

hi <- DE_genes[DE_genes$bnum %in% names,]
