#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("topGO")


library(topGO)
library("org.EcK12.eg.db")

dat = "~/git_repos/HUTI-RNAseq/analysis/DE/2017-06-16-GO-analysis/data/"

fig = "~/git_repos/HUTI-RNAseq/analysis/DE/2017-06-16-GO-analysis/data/"

myDat <- read.csv(file=paste0(dat,"DE-seq-results-all-core-edited.csv"), row.names=1, stringsAsFactors=F)
goDat <- read.csv(file=paste0(dat, "gene_association.ecocyc.csv"), row.names=1, stringsAsFactors=F) # 19361

toRM <- union(which(myDat$sym == "") , which(is.na(myDat$sym)))
myDat1 <- myDat[-toRM,]


iFCok1 <- which( abs(myDat1$log2FoldChange) > 1.5)    # 753 
iFCok2 <- which( abs(myDat1$log2FoldChange) > 0.585)  # 1658 
iFDRok <- which(myDat1$padj < 0.05)           # 1763  
iDE <- intersect(iFDRok,iFCok1)  


geneNam <- as.vector(myDat1$sym) # 2457
DEgeneNam <- as.vector(myDat1$sym[iDE]) # 752

#geneNam1 <- unlist(strsplit(geneNam,","))  # 2776
#DEgeneNam1 <- unlist(strsplit(DEgeneNam,","))  # 118

# Normalize the sets ----
geneNam2 <- geneNam[which(geneNam %in% goDat$Symbol)] # 1804
DEgeneNam2 <- DEgeneNam[which(DEgeneNam %in% goDat$Symbol)]  # 565
goDat2 <- goDat[which(goDat$Symbol %in% geneNam2),]   # 11149

EcoliGenes <- rep(0,length(geneNam2))
EcoliGenes[which(geneNam2 %in% DEgeneNam2  )] <- 1
names(EcoliGenes) <- geneNam2

DEcoli <- names(EcoliGenes)[which(EcoliGenes == 1)]
FacGenes <- as.factor(EcoliGenes)

# make GO to genes lists -----
goBP <- goDat2[which(goDat2$category == "P"),]
goMF <- goDat2[which(goDat2$category == "F"),]
goCC <- goDat2[which(goDat2$category == "C"),]

BPterms <- unique(goBP$GO)  # 1006
MFterms <- unique(goMF$GO)  # 1121
CCterms <- unique(goCC$GO)  # 118

BPlist <- list()#1006
for(ii in 1:length(BPterms)) {
  term <- BPterms[ii]
  genes <- unique(goBP$Symbol[which(goBP$GO == term)])
  BPlist[[ii]] <-  genes
}
names(BPlist) <- BPterms





BPlist2 <- list() #481

for (j in 1:length(BPlist)){
  if (length(BPlist[[j]]) != 1){
    BPlist2 <- append(BPlist2, BPlist[j])
  }
}



MFlist <- list() #1121
for(ii in 1:length(MFterms)) {
  term <- MFterms[ii]
  genes <- unique(goMF$Symbol[which(goMF$GO == term)])
  MFlist[[ii]] <-  genes
}
names(MFlist) <- MFterms

MFlist2 <- list()#355

for (j in 1:length(MFlist)){
  #print(MFlist[[j]])
  if (length(MFlist[[j]]) != 1){
    MFlist2 <- append(MFlist2, MFlist[j])
  }
}

CClist <- list() #118
for(ii in 1:length(CCterms)) {
  term <- CCterms[ii]
  genes <- unique(goCC$Symbol[which(goCC$GO == term)])
  CClist[[ii]] <-  genes
}
names(CClist) <- CCterms

CClist2 <- list()#77

for (j in 1:length(CClist)){
  #print(MFlist[[j]])
  if (length(CClist[[j]]) != 1){
    CClist2 <- append(CClist2, CClist[j])
  }
}


###########

topBP <- new("topGOdata", description = "Ecoli BP", ontology = "BP", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = BPlist2)
topMF <- new("topGOdata", description = "Ecoli MF", ontology = "MF", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = MFlist2)
topCC <- new("topGOdata", description = "Ecoli CC", ontology = "CC", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = CClist2)

BP.genes <- genesInTerm(topBP)
MF.genes <- genesInTerm(topMF)
CC.genes <- genesInTerm(topCC)


# BP -------------------
BP.Fisher.elim <- runTest(topBP, algorithm = "elim", statistic = "fisher")  
BP.Fisher.elim.Table <- GenTable(topBP, elimFisher=BP.Fisher.elim, topNodes=1000)
BP.Fisher.elim.Table <- BP.Fisher.elim.Table[which(BP.Fisher.elim.Table$elimFisher < 0.05),]  # 70

if(nrow(BP.Fisher.elim.Table) > 0) { 
  goIDs <- BP.Fisher.elim.Table$GO.ID
  genesInTerm <- sapply(goIDs, function(x) paste(BP.genes[[which(names(BP.genes) == x)]], collapse=","))
  DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(BP.genes[[which(names(BP.genes) == x)]], DEgeneNam1), collapse=","))
  aDF <- data.frame(genesInTerm,DEGenesInTerm)
  bDF <- cbind(BP.Fisher.elim.Table,aDF)
  write.csv(bDF, file=paste0(dat, "BP.Fisher.elim.Table.csv") )
}


          
          
# MF -------------------
MF.Fisher.elim <- runTest(topMF, algorithm = "elim", statistic = "fisher")  
MF.Fisher.elim.Table <- GenTable(topMF, elimFisher=MF.Fisher.elim, topNodes=651)
MF.Fisher.elim.Table <- MF.Fisher.elim.Table[which(MF.Fisher.elim.Table$elimFisher < 0.05),]  #23


if(nrow(MF.Fisher.elim.Table) > 0) { 
  goIDs <- MF.Fisher.elim.Table$GO.ID
  genesInTerm <- sapply(goIDs, function(x) paste(MF.genes[[which(names(MF.genes) == x)]], collapse=","))
  DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(MF.genes[[which(names(MF.genes) == x)]], DEgeneNam1), collapse=","))
  aDF <- data.frame(genesInTerm,DEGenesInTerm)
  bDF <- cbind(MF.Fisher.elim.Table,aDF)
  write.csv(bDF, file=paste0(dat, "MF.Fisher.elim.Table.csv") )
}

# CC -------------------

CC.Fisher.elim <- runTest(topCC, algorithm = "elim", statistic = "fisher")  
CC.Fisher.elim.Table <- GenTable(topCC, elimFisher=CC.Fisher.elim, topNodes=150)
CC.Fisher.elim.Table <- CC.Fisher.elim.Table[which(CC.Fisher.elim.Table$elimFisher < 0.05),]  #9
#if(nrow(CC.Fisher.elim.Table) > 0) { write.csv(CC.Fisher.elim.Table, file="CC.Fisher.elim.Table.csv") }
if(nrow(CC.Fisher.elim.Table) > 0) { 
  goIDs <- CC.Fisher.elim.Table$GO.ID
  genesInTerm <- sapply(goIDs, function(x) paste(CC.genes[[which(names(CC.genes) == x)]], collapse=","))
  DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(CC.genes[[which(names(CC.genes) == x)]], DEgeneNam1), collapse=","))
  aDF <- data.frame(genesInTerm,DEGenesInTerm)
  bDF <- cbind(CC.Fisher.elim.Table,aDF)
  write.csv(bDF, file=paste0(dat,"CC.Fisher.elim.Table.csv") )
}
# -----







