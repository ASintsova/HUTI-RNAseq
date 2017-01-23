#Dependencies

library(DESeq2)
library (pheatmap)
library(RColorBrewer)
library (ggplot2)

#RNAseq Counts Analysis

#INPUT: experiment info + count data

RNAseq_analysis <- function( count_data, experiment_info, design_par){

        count_dds  <- DESeqDataSetFromMatrix(countData = as.matrix(count_data), colData = experiment_info,
                                            design = design_par)
        return (count_dds)
}



# draw heatmap and PCA analysis using DESeq2, takes in DESeq object

#NEEDS TO BE TROUBLESHOOTED!!!!!
RNAseq_cluster_DESeq <- function(count_data_dds, design_par){
        dds_transformed <- rlog(count_data_dds, blind = FALSE) #transform
        core_dists <- dist(t(assay(dds_transformed)))# calculated distances
        dists_matrix <- as.matrix(core_dists)
        par = substring(design_par, 1)[2]
        #par = design_par
        print (par)
        rownames(dists_matrix)<- paste(dds_transformed$STRAIN, dds_transformed$par, sep = "_")
        colnames(dists_matrix) <- NULL
        
        colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
        pheatmap(dists_matrix, #draw heatmap
                 clustering_distance_rows=core_dists,
                 clustering_distance_cols=core_dists,
                 col=colors)
        plotPCA(dds_transformed, intgroup = c(par)) #draw PCA
        
}
