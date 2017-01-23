#Dependencies

library(DESeq2)


#RNAseq Counts Analysis

#INPUT: experiment info + count data

RNAseq_analysis <- function(experiment_info, count_data, design_par = colnames(experiment_info)[1]){
        count_data <- as.matrix(count_data)
        (count_data_dds <- DESeqDataSetFromMatrix(countData = count_data, colData = experiment_info,
                                            design = ~design_par))

        
        
}
# draw heatmap and PCA analysis using DESeq2, takes in DESeq object
#NEEDS TO BE TROUBLESHOOTED
RNAseq_exploratory_analysis_and_visualization <- function(count_data_dds, design_par){
        dds_transformed <- rlog(count_data_dds, blind = FALSE) #transform
        core_dists <- dist(t(assay(dds_transformed)))# calculated distances
        dists_matrix <- as.matrix(core_dists)
        rownames(dists_matrix)<- paste(dds_transformed$STRAIN, dds_transformed$design_par, sep = "_")
        colnames(dists_matrix) <- NULL
        
        colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
        pheatmap(dists_matrix, #draw heatmap
                 clustering_distance_rows=core_dists,
                 clustering_distance_cols=core_dists,
                 col=colors)
        plotPCA(core_dds_transformed, intgroup = c(design_par)) #draw PCA
        
}
