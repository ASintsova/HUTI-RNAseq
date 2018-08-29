import pandas as pd
import os

# todo remove this!!!
genome = "HM01"
treatment = "UR"
count_path = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/data/HM01_UR_trimmed_sorted_counts"
gff_path = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes/data/gff_files/HM01.gff"
flagstat_summary = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/counts/flagstat/flagstat_summary.txt"
crossRef = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes/data/genomes_homologues/intersect_pan_CM_t0/pangenome_matrix_t0_crossRef.csv"

def addGeneLength(genome, treatment, count_path, gff_path):

    # Read in counts
    co = pd.read_csv(count_path,sep ="\t", index_col=0, header=None)
    co.columns = ["{}_{}_counts".format(genome, treatment)]
    to_drop = [ind for ind in list(co.index) if "__" in ind]
    co.drop(to_drop, inplace=True)
    # Add gene length to counts
    gene_length = {}
    with open (gff_path, "r") as gff:
        for line in gff:
            if "CDS" not in line:
                continue
            try:
                gene = line.split(";")[0].split("gene_id=")[1]
                gene_length[gene] = abs(int(line.split()[4])-int(line.split()[3]))
            except:
                pass
    co["gene_length"] = pd.Series(gene_length)
    return co



def addRPKM(genome, treatment, counts_with_gene_length, config):
    flagstat = config.get("FLAG", "path")
    flag = config.get("FLAG", "format").replace("genome", genome).replace("sample", treatment)
    with open(flagstat, "r") as fs:
        for line in fs:
            if flag in line:
                N = int(line.split(',')[4])
    counts_with_gene_length["{}_{}_RPKM".format(genome, treatment)] = round((10**9*counts_with_gene_length["{}_{}_counts".format(genome, treatment)])/(N*counts_with_gene_length["gene_length"]),2)
    counts_with_gene_length.drop(["gene_length"], axis=1, inplace=True)
    return counts_with_gene_length


print(os.path.isfile(count_path))
print(os.path.isfile(gff_path))
co = addGeneLength(genome, treatment, count_path, gff_path)

def joinAll(crossRef, counts_with_gene_length):

    ortho_df = pd.read_csv(crossRef, index_col=0)
    print(ortho_df.head())
    print(co.head())
    final = ortho_df.join(co, on="HM01")
    print(final.dropna().head())

joinAll(crossRef, co)

