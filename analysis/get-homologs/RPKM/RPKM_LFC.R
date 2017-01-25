#RPKM analysis

# Want to transform
# look at means and FC

rpkms <- read.csv("~/git_repos/HUTI-RNAseq/analysis/get-homologs/RPKM/combined_RPKM.csv", header = TRUE, row.names = 1)
rpkms_log <- log2(rpkms+1)
rpkms_means <- rowMeans(rpkms_log, na.rm = TRUE)
odd <- seq(1,26, 2)
even <- seq(2,26,2)
urine_RPKM <- rpkms_log[,odd]
uti_RPKM <- rpkms_log[,even]
ur_rpkms_means <- rowMeans(rpkms_log[,odd], na.rm = TRUE)
uti_rpkms_means <- rowMeans(rpkms_log[,even], na.rm = TRUE)


rpkms_means_LFC <- cbind(rpkms_log, ur_rpkms_means, uti_rpkms_means)

rpkms_LFC <- subset(rpkms_means_LFC, ur_rpkms_means!=0 & uti_rpkms_means!=0)

rpkms_LFC <- transform(rpkms_LFC, LFC = uti_rpkms_means/ur_rpkms_means)


rpkms_LFC <- cbind(rpkms_LFC, na_num = rowSums(is.na(rpkms_LFC)))

#x <- rpkms_means_LFC[order(rpkms_means_LFC$na_num, -rpkms_means_LFC$LFC),]

y <- rpkms_LFC[order(-rpkms_LFC$LFC),]

#write.csv(an, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/RPKM/ur_vs_uti_means.csv", row.names=TRUE)
write.csv(y, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/RPKM/ur_vs_uti_means.csv", row.names=TRUE)

r.means <- read.csv("~/git_repos/HUTI-RNAseq/analysis/get-homologs/RPKM/2017-01-25-ur_vs_uti_RPKM_means-edited.csv", header = TRUE, row.names = 1)


r.means2 <- subset(r.means, (LFC > 1| LFC <0.5) & na_num > 0)



rpkm.means <- as.data.frame(r.means2[, c(29)])
row.names(rpkm.means) <- make.names(r.means2[, 31], TRUE)

rpkm.means <- apply(rpkm.means, 2, function(x) ifelse(x <1, -1/x, x))

symbols <- row.names(rpkm.means)

head(symbols)


ec.eg <- id2eg(ids = symbols, category = gene.idtype.list[1], org = "EcK12")

colnames(ec.eg) <- c("gene_name", "ENSMBL")
head(ec.eg)

gene_list <- cbind(ec.eg[,2], as.numeric(rpkm.means[,1])) # includes ensmbl number thing
colnames(gene_list) <- c("ENSMBL", colnames(rpkm.means))
head(gene_list)

tail(gene_list)


gene_list <- gene_list[complete.cases(gene_list),]

analyze_list <- as.matrix(as.numeric(gene_list[,2]))
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


