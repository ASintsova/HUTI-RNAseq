csvfile <- file.path("GTFs", "sample_info.csv")
(sample_info <- read.csv(csvfile,row.names=1))


strain_nums <- c(1, 3, 6, 7, 14, 17, 43, 54, 56, 57, 66, 68, 86)
#strain_nums <- c(1)
res <- c()
for (i in strain_nums){
    re <- countReads(i)
    re
    write.csv(re, paste0("HM", i, "_gene_counts"))
    
    
}


countReads <- function (num)
{
    library("Rsamtools")
    library("GenomicFeatures")
    library("GenomicAlignments")
    library("BiocParallel")
    filenames <- file.path("GTFs", paste0(levels(sample_info$Condition), num, "_SE_aln_sort.bam"))
    file.exists(filenames)
    bamfiles <- BamFileList(filenames, yieldSize=2000000)
    gtffile <- file.path('core_genes_by_genome', paste0("HM", num, "_genes.df"))
    file.exists(gtffile)
    gtf <- read.table(gtffile, header = TRUE)
    ebg <- makeGRangesFromDataFrame(gtf, keep.extra.columns = TRUE, seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand")
    
    se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                            mode="Union",
                            singleEnd=TRUE,
                            ignore.strand=FALSE )
    
    rownames(se)<-ebg$feature
    assay(se)
   
}


