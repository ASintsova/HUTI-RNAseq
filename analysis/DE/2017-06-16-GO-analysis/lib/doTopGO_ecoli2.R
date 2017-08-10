
library("topGO")
library("org.EcK12.eg.db")

# get the gene list and convert to entrez
myDat <- read.csv(file="/Users/rtagett/UMCore/Clients/ZNotCurrent/Mobley/ForComps/SPARTA_results/DE/TosR-Emp_expression_noEBG_Annot.csv", row.names=1, stringsAsFactors=F)
goDat <- read.csv(file="/Users/rtagett/UMCore/Clients/ZNotCurrent/Mobley/ForComps/FunctionalAnalysis/GO/gene_association.ecocyc.csv", row.names=1, stringsAsFactors=F) # 19361

toRM <- union(which(myDat$UpdatedSymbol == "") , which(is.na(myDat$UpdatedSymbol)))
myDat1 <- myDat[-toRM,]


iFCok1 <- which( abs(myDat1$logFC) > 1.5)    # 151 (83) - after removing those with no gene Symbol
iFCok2 <- which( abs(myDat1$logFC) > 0.585)  # 842 (606)
iFDRok <- which(myDat1$FDR < 0.05)           # 200 (113) 
iDE <- intersect(iFDRok,iFCok2)  # 200 (113)

geneNam <- as.vector(myDat1$UpdatedSymbol) # 2538
DEgeneNam <- as.vector(myDat1$UpdatedSymbol[iDE]) # 113

geneNam1 <- unlist(strsplit(geneNam,","))  # 2776
DEgeneNam1 <- unlist(strsplit(DEgeneNam,","))  # 118

# Normalize the sets ----
geneNam2 <- geneNam1[which(geneNam1 %in% goDat$Symbol)] # 1987
DEgeneNam2 <- DEgeneNam1[which(DEgeneNam1 %in% goDat$Symbol)]  # 70
goDat2 <- goDat[which(goDat$Symbol %in% geneNam2),]   # 13378

EcoliGenes <- rep(0,length(geneNam2))
EcoliGenes[which(geneNam2 %in% DEgeneNam2  )] <- 1
names(EcoliGenes) <- geneNam2

DEcoli <- names(EcoliGenes)[which(EcoliGenes == 1)]
FacGenes <- as.factor(EcoliGenes)

# make GO to genes lists -----
goBP <- goDat2[which(goDat2$category == "P"),]
goMF <- goDat2[which(goDat2$category == "F"),]
goCC <- goDat2[which(goDat2$category == "C"),]

BPterms <- unique(goBP$GO)  # 1110
MFterms <- unique(goMF$GO)  # 1228
CCterms <- unique(goCC$GO)  # 122

BPlist <- list()
for(ii in 1:length(BPterms)) {
    term <- BPterms[ii]
    genes <- unique(goBP$Symbol[which(goBP$GO == term)])
    BPlist[[ii]] <-  genes
}
names(BPlist) <- BPterms

MFlist <- list()
for(ii in 1:length(MFterms)) {
    term <- MFterms[ii]
    genes <- unique(goMF$Symbol[which(goMF$GO == term)])
    MFlist[[ii]] <-  genes
}
names(MFlist) <- MFterms

CClist <- list()
for(ii in 1:length(CCterms)) {
    term <- CCterms[ii]
    genes <- unique(goCC$Symbol[which(goCC$GO == term)])
    CClist[[ii]] <-  genes
}
names(CClist) <- CCterms


###########

topBP <- new("topGOdata", description = "Ecoli BP", ontology = "BP", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = BPlist)
topMF <- new("topGOdata", description = "Ecoli MF", ontology = "MF", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = MFlist)
topCC <- new("topGOdata", description = "Ecoli CC", ontology = "CC", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = CClist)

BP.genes <- genesInTerm(topBP)
MF.genes <- genesInTerm(topMF)
CC.genes <- genesInTerm(topCC)


# BP -------------------
BP.Fisher.elim <- runTest(topBP, algorithm = "elim", statistic = "fisher")  
BP.Fisher.elim.Table <- GenTable(topBP, elimFisher=BP.Fisher.elim, topNodes=1000)
BP.Fisher.elim.Table <- BP.Fisher.elim.Table[which(BP.Fisher.elim.Table$elimFisher < 0.05),]  # 37

if(nrow(BP.Fisher.elim.Table) > 0) { 
    goIDs <- BP.Fisher.elim.Table$GO.ID
    genesInTerm <- sapply(goIDs, function(x) paste(BP.genes[[which(names(BP.genes) == x)]], collapse=","))
    DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(BP.genes[[which(names(BP.genes) == x)]], DEgeneNam1), collapse=","))
    aDF <- data.frame(genesInTerm,DEGenesInTerm)
    bDF <- cbind(BP.Fisher.elim.Table,aDF)
    write.csv(bDF, file="BP.Fisher.elim.Table.csv") 
}


# MF -------------------
MF.Fisher.elim <- runTest(topMF, algorithm = "elim", statistic = "fisher")  
MF.Fisher.elim.Table <- GenTable(topMF, elimFisher=MF.Fisher.elim, topNodes=1000)
MF.Fisher.elim.Table <- MF.Fisher.elim.Table[which(MF.Fisher.elim.Table$elimFisher < 0.05),]  # 46

if(nrow(MF.Fisher.elim.Table) > 0) { 
    goIDs <- MF.Fisher.elim.Table$GO.ID
    genesInTerm <- sapply(goIDs, function(x) paste(MF.genes[[which(names(MF.genes) == x)]], collapse=","))
    DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(MF.genes[[which(names(MF.genes) == x)]], DEgeneNam1), collapse=","))
    aDF <- data.frame(genesInTerm,DEGenesInTerm)
    bDF <- cbind(MF.Fisher.elim.Table,aDF)
    write.csv(bDF, file="MF.Fisher.elim.Table.csv") 
}

# CC -------------------

CC.Fisher.elim <- runTest(topCC, algorithm = "elim", statistic = "fisher")  
CC.Fisher.elim.Table <- GenTable(topCC, elimFisher=CC.Fisher.elim, topNodes=200)
CC.Fisher.elim.Table <- CC.Fisher.elim.Table[which(CC.Fisher.elim.Table$elimFisher < 0.05),]  # 2
if(nrow(CC.Fisher.elim.Table) > 0) { write.csv(CC.Fisher.elim.Table, file="CC.Fisher.elim.Table.csv") }
if(nrow(CC.Fisher.elim.Table) > 0) { 
    goIDs <- CC.Fisher.elim.Table$GO.ID
    genesInTerm <- sapply(goIDs, function(x) paste(CC.genes[[which(names(CC.genes) == x)]], collapse=","))
    DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(CC.genes[[which(names(CC.genes) == x)]], DEgeneNam1), collapse=","))
    aDF <- data.frame(genesInTerm,DEGenesInTerm)
    bDF <- cbind(CC.Fisher.elim.Table,aDF)
    write.csv(bDF, file="CC.Fisher.elim.Table.csv") 
}
# -----







