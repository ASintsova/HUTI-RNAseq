

#calculating K statistic, saving values of K-statistic P-value and normalized counts
#returns a list of 2 dataframes 

library(phytools)
library(DESeq2)



calculateKForDataFrame <- function(trait_matrix, tree, filename) {
        print ("K stat and P values can be accessed via $phylosig_K, matrix values: $phylosig_gene_counts")
        K_values <- c()
        P_values <- c()
        genes_with_signal <- c()
        mat <- matrix(NA, nrow=nrow(trait_matrix), ncol=3)
        
        for (i in 1:nrow(trait_matrix)){
                
                gene_name = rownames(trait_matrix[i,])
                expression_vector <- as.numeric(trait_matrix[i,])
                names(expression_vector) <- c("HM1", "HM3", "HM6", "HM7", "HM14", "HM17", "HM43", "HM54", "HM56", "HM57", "HM66", "HM68", "HM86")
                K <- phylosig(tree, expression_vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                              control=list())$K
                P <- phylosig(tree, expression_vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                              control=list())$P
                if (P <= 0.01){
                        genes_with_signal <- append(genes_with_signal, c(gene_name, K, P))
                        mat[i,]<- c(gene_name, K, P)
                }
                K_values <- append(K_values, K)
                P_values <- append(P_values, P)
        }
        
        mat <- mat[complete.cases(mat),]
        
        phylo_sig <- mat[,c(2,3)]
        rownames(phylo_sig) <- mat[,1]
        signal_names <- mat[,1]
        colnames(phylo_sig)<- c("K-statistic", "P-value")
        phylosig_traits <- trait_matrix[signal_names,]
        
        phy <- list("phylosig_K" = as.data.frame(phylo_sig), "phylosig_gene_counts" = phylosig_traits)
        fn = paste0(DATAPATH, filename)
        write.csv(phylo_sig, fn, row.names=TRUE)
        return (phy)
}
