dat = "~/git_repos/HUTI-RNAseq/analysis/DE/2017-06-19-exploration/data/"

fig = "~/git_repos/HUTI-RNAseq/analysis/DE/2017-06-19-exploration/figures/"

library(pheatmap)

myDat = read.csv(file=paste0(dat,"DE-seq-results-all-core-edited.csv"), row.names=1, stringsAsFactors=F)

info <- read.csv(file=paste0(dat,"Patient_meta_info.csv"), row.names=1)

topup <- head(myDat[order(myDat$log2FoldChange, decreasing = TRUE),], 25)
write.csv(topup, paste0(dat, "top-upregulated-genes.csv"))


#

topdown<- head(myDat[order(myDat$log2FoldChange),], 25)


## Attempted to use topGO to perform enrichment analysis on regulons -> that was not successful



counts = read.csv(file=paste0(dat,"normalized_counts_core.csv"), row.names=1, stringsAsFactors=F)

counts_sym <- merge(counts, myDat, by= "row.names", all.x = TRUE)

counts2 <- counts_sym[,c(2:27,29,34)]

counts3 <- cbind(counts2[,seq(1,26,2)], counts2[,c(seq(2,26,2),28)])
rownames(counts3)<-rownames(counts2)

#Gluconate

gluconate <- c("gntU", "gntK", "gntR")
annot <- info[,c(2,5)]
gluc_counts <- counts2[which(counts2$sym %in% gluconate),]
names <- gluc_counts$sym
gluc_counts <- log2(gluc_counts[,1:26] +1)
rownames(gluc_counts)<- names
pheatmap(gluc_counts, annotation_col = annot,cellheight = 15, cellwidth = 15,
         cluster_rows = FALSE)


#Biotin

biotin <- c("bioB", "bioC", "bioD","bioF")
biotin_counts <- counts3[which(counts3$sym %in% biotin),]
names <- biotin_counts$sym
biotin_counts <- log2(biotin_counts[,1:26] +1)
rownames(biotin_counts)<- names
pheatmap(biotin_counts, annotation_col = annot,cellheight = 15, cellwidth = 15,
         cluster_rows = FALSE, cluster_cols = FALSE)


#making a scatter plot


biotinData <- list()
biotinData$id <- rep(rownames(t(biotin_counts)), 4)
biotinData$strain <- rep(c("HM1", "HM3", "HM6", "HM7", "HM14", "HM17", "HM43",
                           "HM54", "HM56", "HM57", "HM66", "HM68", "HM86"), 4)
biotinData$media <- rep(c(rep("UR", 13), rep("UTI", 13)),4)
biotinData$gene <- c(rep("bioB", 26), rep("bioF", 26), rep("bioC", 26), rep("bioD", 26))
biotinData$counts<- c(as.numeric(biotin_counts[1, ]), as.numeric(biotin_counts[2,]),as.numeric(biotin_counts[3,]),
                      as.numeric(biotin_counts[4,]))
biotinData <- as.data.frame(biotinData)
                  
library(ggplot2)
ggplot(biotinData, aes(x = media, y = counts, color = strain))  +
  geom_jitter(position=position_jitter(width=.1, height=0))+
  geom_line(aes(group = strain),
            alpha = 0.5, colour = "darkgrey")+
  facet_wrap(~ gene, ncol = 4)



#glycogen

glycogen <- c("glgA", "glgC", "glgP")
glycogen_counts <- counts3[which(counts3$sym %in% glycogen),]
names <- glycogen_counts$sym
glycogen_counts <- log2(glycogen_counts[,1:26] +1)
rownames(glycogen_counts)<- names
pheatmap(glycogen_counts, annotation_col = annot,cellheight = 15, cellwidth = 15,
         cluster_rows = FALSE, cluster_cols = FALSE)
#curli



curli <- c("csgA", "csgB", "csgF", "csgG", "csgE", "bssR", 
           "glgS", "rpoS", "bolA")
curli_counts <- counts3[which(counts3$sym %in% curli),]
names <- curli_counts$sym
curli_counts <- log2(curli_counts[,1:26] +1)
rownames(curli_counts)<- names
pheatmap(curli_counts, annotation_col = annot,cellheight = 15, cellwidth = 15,
         cluster_rows = FALSE, cluster_cols = FALSE)




#wzb
curli <- c("wzb", "wza")
curli_counts <- counts3[which(counts3$sym %in% curli),]
names <- curli_counts$sym
curli_counts <- log2(curli_counts[,1:26] +1)
rownames(curli_counts)<- names
pheatmap(curli_counts, annotation_col = annot,cellheight = 15, cellwidth = 15,
         cluster_rows = FALSE)


#stress response

stress <- c("rpoS", "hdeA", "cspD", "yhcJ", "tnaB", "yjiY")
stress_counts <- counts3[which(counts3$sym %in% stress),]
names <- stress_counts$sym
stress_counts <- log2(stress_counts[,1:26] +1)
rownames(stress_counts)<- names
pheatmap(stress_counts, annotation_col = annot,cellheight = 15, cellwidth = 15,
         cluster_rows = FALSE)

