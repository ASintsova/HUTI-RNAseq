from pathlib import Path
import rpy2
import utilities as ut
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

config_dict = ut.process_config()




def get_raw_counts_for_R(config_dict):

    counts_dir = config_dict['data']['counts']
    orth_matrix = pd.read_csv(config_dict['data']['ortholog_table'], index_col=0)
    # Gather raw counts for R
    dfl = []
    for f in Path(counts_dir).iterdir():
        df = pd.read_table(f, index_col=0, header=None, names=[f.stem.split("_counts")[0]])
        dfl.append(df)
    raws = pd.concat(dfl, axis=1, sort=False)
    raw_counts = ut.get_tpms_for_prokkas(orth_matrix, raws)
    raw_counts = raw_counts[sorted(raw_counts.columns)]
    raw_counts.to_csv('shared_data/raw_core_counts.csv')



deseq2 = importr('DESeq2')


# Defining
DEseqFromMat = robjects.r['DESeqDataSetFromMatrix']
DESeq = robjects.r['DESeq']
estimateSizeFactors = robjects.r['estimateSizeFactors']
get_norm_counts = robjects.r['counts']
results = robjects.r['results']
formula = robjects.r['as.formula']
hd = robjects.r['head']
pData = robjects.r['pData']
DF = robjects.r['data.frame']

robjects.r('''
        process_results <- function (res, prefix, outdir){
        resSig <- subset(res, (padj < 0.05)& (abs(log2FoldChange) > 2))
        write.csv(resSig[order(resSig$log2FoldChange, decreasing = TRUE),], 
                  file.path(outdir, paste0(prefix, "_significant.csv")),
                  row.names = TRUE, quote = FALSE)
        write.csv(res[order(res$log2FoldChange, decreasing = TRUE),], 
        file.path(outdir, paste0(prefix, "_all.csv")),row.names = TRUE, quote = FALSE)
        }
        ''')

process_results = robjects.r['process_results']


def get_deseq_norm_counts(counts, colData, variable):

    dds = DEseqFromMat(countData=counts, colData=colData, design=formula('~ {}'.format(variable)))
    analysis = DESeq(dds)
    esf = estimateSizeFactors(analysis)
    norm_counts = get_norm_counts(esf, normalized=True)
    df = pd.DataFrame(norm_counts, index=counts.index, columns=counts.columns)
    return df



def de_analysis(counts, colData, variable, treat, untreat, prefix, lfcThr=0):
    dds = DEseqFromMat(countData=counts, colData=colData,
                       design=formula('~ {}'.format(variable)))

    analysis = DESeq(dds)

    v = rpy2.robjects.vectors.StrVector([variable, treat, untreat])

    res = results(analysis, contrast=v, lfcThreshold=lfcThr)
    df = pd.DataFrame(DF(res))
    return df

def process_results(res, prefix, out_dir):
    pass



if __name__ == '__main__':
    counts_file = config_dict['data']['raw_core_counts']
    colData_file = config_dict['data']['patient_info_r']
    counts = pd.read_csv(counts_file, index_col=0)
    print(counts.head())
    colData = pd.read_csv(colData_file, index_col=0)
    df = de_analysis(counts, colData, "MEDIA", "PATIENT", "URINE", "TEST", lfcThr=0)
    print(df.head())












