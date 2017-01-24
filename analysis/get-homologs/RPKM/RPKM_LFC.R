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
ur_rpkms_means[uti_rpkms_means == 0] <- NA

rpkms_means_LFC <- cbind(rpkms_log, ur_rpkms_means, uti_rpkms_means)

rpkms_means_LFC <- transform(rpkms_means_LFC, LFC = uti_rpkms_means/ur_rpkms_means)


rpkms_means_LFC <- cbind(rpkms_means_LFC, na_num = rowSums(is.na(rpkms_log)))

x <- rpkms_means_LFC[order(rpkms_means_LFC$na_num, -rpkms_means_LFC$LFC),]

write.csv(an, "~/git_repos/HUTI-RNAseq/analysis/get-homologs/RPKM/ur_vs_uti_means.csv", row.names=TRUE)

