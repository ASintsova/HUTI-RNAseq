library(phytools)
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/KStatFromTraitMatrix.R")
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LambdaFromTraitMatrix.R")
source ("~/git_repos/HUTI-RNAseq/analysis/evolutionaryModels/lib/LoadHTseqCountsandTree.R")


all_data <- LoadHTSeqCountsandTree()


tree <- all_data$tree

ace <- as.numeric(all_data$uti_counts_l["17669_aceE",])
names(ace)<- c("HM1", "HM3", "HM6", "HM7", "HM14", "HM17", "HM43", "HM54", "HM56", "HM57", "HM66", "HM68", "HM86")
phenogram(tree, ace, spread.labels = TRUE)

obj <- contMap(tree, ace)


ace_ur <- as.numeric(all_data$urine_counts_l["17669_aceE",])
names(ace_ur)<- c("HM1", "HM3", "HM6", "HM7", "HM14", "HM17", "HM43", "HM54", "HM56", "HM57", "HM66", "HM68", "HM86")

obj2 <- contMap(tree, ace_ur)
