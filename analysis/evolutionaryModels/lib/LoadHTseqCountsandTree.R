#Load in the data for phylosignal analysis


LoadHTSeqCountsandTree <- function(counts_filename="~/git_repos/HUTI-RNAseq/analysis/DE/2017-01-30-htseq-counts/normalized_counts_core.csv",
                                   tree_filename = NULL, 
                                   metadata = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/data/Patient_meta_info.csv"){
        odd <- seq(1,26,2)
        even <- seq(2, 26, 2)
        
        experiment_info <- read.csv(metadata, header = TRUE, row.names=1)
        
        counts <- read.csv(counts_filename, header = TRUE, row.names = 1)
        if (is.null(tree_filename)){
                full_tree <- read.tree(text ="((HM17:0.00428,HM57:0.00876)1.000:0.08988,(HM65:0.00824,HM68:0.00963)1.000:0.09530,(((HM54:0.08292,HM66:0.07434)1.000:0.04338,(HM43:0.11275,((HM60:0.24376,((HM7:0.07135,(HM14:0.08006,HM3:0.07130)1.000:0.02392)1.000:0.40723,(HM46:0.23512,(HM26:0.19620,(HM1:0.01633,HM69:0.01337)1.000:0.17697)1.000:0.04329)1.000:0.10117)1.000:0.07596)1.000:0.33315,(HM56:0.08289,HM6:0.08508)1.000:0.08846)1.000:0.04832)1.000:0.01448)1.000:0.01347,(HM86:0.10068,HM27:0.11440)0.970:0.01392)1.000:0.01738);")
                #focusing on 13 strains 
                tree <- drop.tip(full_tree, c("HM65", "HM27", "HM26", "HM69", "HM46", "HM60"))
        }
        else{
                tree <- read.tree(tree_filename)
        }


        urine_counts <- counts[, odd]
        uti_counts <- counts[,even]
        urine_counts_l <- log2(urine_counts +1)
        uti_counts_l <- log2(uti_counts +1)
        urine_info <- experiment_info[odd,]
        uti_info <- experiment_info[even,]
        data <- list("counts" = counts, "urine_counts" = urine_counts, "uti_counts" = uti_counts, "urine_counts_l" = urine_counts_l, "uti_counts_l" = uti_counts_l,
                     "urine_info" = urine_info, "uti_info" = uti_info, "tree" = tree, "exp_info" = experiment_info)
        
        return (data)         
}

d <- LoadHTSeqCountsandTree()

