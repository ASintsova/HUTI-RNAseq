
require(graph)
require(ROntoTools)


myDat <- read.csv(file="/Users/rtagett/UMCore/Clients/Mobley/Sparta_Results/DE/TosR-Emp_expression.csv", row.names=1, stringsAsFactors=F)

ref <- paste("ecc:",myDat$geneID, sep="")
myDat <- cbind(ref,myDat)
# kpg_ecc.RData

#--------------


#kpg <- keggPathwayGraphs("ecc",nodeOnlyGraphs=TRUE, updateCache = TRUE,  verbose = TRUE, relPercThresh=0)  # expect 118 pathways
#kpgBAK <- kpg
#save(kpgBAK, file="kpg_ecc.RData")

load(file="kpg_ecc.RData")
kpg <- kpgBAK

# kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
#  edgeWeightByType = list(activation = 1, inhibition = -1,
#  expression = 1, repression = -1),
#  defaultWeight = 0)

 kpn <- keggPathwayNames("ecc")



iOK <- intersect(which(abs(myDat$logFC) > 0.585) , which(myDat$FDR < 0.05 ) )
fc <- myDat$logFC[iOK]
names(fc) <- ref[iOK]
pv <- myDat$FDR[iOK]
names(pv) <- ref[iOK]

    kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)

    peRes <- pe(x = fc, graphs = kpg, ref = ref, nboot = 200, verbose = FALSE)

    res <- Summary(peRes, pathNames = kpn, totalAcc = TRUE, totalPert = TRUE, normalize = TRUE,pPert = TRUE, 
		pAcc = TRUE, pORA = TRUE, comb.pv = c("pPert", "pORA"), 
		comb.pv.func = compute.fisher,order.by = "pComb", adjust.method = "fdr")

    res <- res[which(res$pComb.fdr > 0), ]


    upGenes <- c() ; downGenes <- c() ; nbGenesKegg <- c()
    for(jj in 1:nrow(res)) {
	path <- as.vector(res$pathNames[jj])
	pw <- rownames(res)[jj]
	p <- peRes@pathways[[pw]]
	g <- layoutGraph(p@map, layoutType = "dot")
	gDF <- myDat[which(myDat$ref %in% nodes(g)),]
	up <- paste(sort(as.vector(gDF$geneID[which(gDF$logFC > 0)])),collapse=",")
        down <- paste(sort(as.vector(gDF$geneID[which(gDF$logFC < 0)])),collapse=",")
	upGenes <- append(upGenes , up)
	downGenes <- append(downGenes, down)
	nbGenesKegg <- append(nbGenesKegg, length(nodes(g)) )
    }

    zDF <- data.frame(upGenes,downGenes,nbGenesKegg)
    res <- cbind(res,zDF)
    
    write.csv(res,file="TosR-Emp_KeggResults.csv")



#  poo <- tryCatch.W.E(peNodeRenderInfo(p))

# p <- peRes@pathways[["path:hsa05216"]]
# g <- layoutGraph(p@map, layoutType = "dot")
# graphRenderInfo(g) <- list(fixedsize = FALSE)
# edgeRenderInfo(g) <- peEdgeRenderInfo(p)
# nodeRenderInfo(g) <- peNodeRenderInfo(p)
# renderGraph(g)








