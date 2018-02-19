"""
    Given a clinical genome (ex. HM01), reference(s), condition (UR/UTI), seq.info (-/seq1/seq2) and get_hom output
    provide counts, RPKMs, and names of orthologs in the reference
    Config file has to have paths to gff files, get_homologues output, htseq-count output directory

"""
import os
import configparser
import pandas as pd
from modules import generating_OGC_crossref as cr
from modules import parsing_gbk as pg
import datetime as datetime

# 1. Check if CrossRef has been generated


#config = configparser.ConfigParser()
#config.read("lib/config")


def get_gene_id(genome, sample, ref_genomes, config):
    GETHOMS = config.get("GETHOMS", "path") # dir
    HTSEQ = config.get("HTSEQ", "path") # dir
    FLAG = config.get("FLAG", "path") # file

    if os.path.isfile(config.get("GFF", "path")):
        gff_path = config.get("GFF", "path")
    else:
        GFF = config.get("GFF", "path") # dir
        gff_file = config.get("GFF", "format").replace("genome", genome)
        gff_path = os.path.join(GFF, gff_file)
    print("Checking get_homologues output...")
    if not checkGETHOMS(genome, ref_genomes, GETHOMS):
        return "get_homologues output does not include all of these genomes, exiting"
    else:
        path = os.path.join(GETHOMS, "intersect_pan{}".format(config.get("GETHOMS", "suffix")))
        assert os.path.isdir(path)
        matrix = os.path.join(path, "pangenome_matrix{}.tab".format(config.get("GETHOMS", "matrix_suffix")))
        assert os.path.isfile(matrix)
        if isCrossRef(matrix):
            pgm = pd.read_csv((matrix.split(".")[0] + "_crossRef.csv"), index_col=0)
        else:
            print("No crossRef file, generating")
            pgm_file = cr.crossRefFromPanGenome(matrix, path)
            gbk, genome = cft073()
            pg.replaceProteinIdWithLocusTags(gbk, genome, pgm_file)
            gbk, genome = MG1655()
            pg.replaceProteinIdWithLocusTags(gbk, genome, pgm_file)
            pgm = pd.read_csv(pgm_file) #only need one argument really
        # Get orthologs
        ortho_df = get_orthologs(pgm, genome, ref_genomes)

        # Get counts
        count_file = config.get("HTSEQ", "format").replace("genome", genome).replace("sample", sample)
        print(count_file)
        count_path = os.path.join(HTSEQ, count_file)
        print(os.path.isfile(count_path))
        #gff_file = config.get("GFF", "format").replace("genome", genome)
        #gff_path = os.path.join(GFF, gff_file)

        # Add counts
        count_df = addCounts(genome, sample, ortho_df, count_path, gff_path)

        # Get total # of reads
        flag_sum = os.path.join(FLAG, "flagstat_summary.txt")
        flag = config.get("FLAG", "format").replace("genome", genome).replace("sample", sample)
        with open(flag_sum, "r") as fs:
            for line in fs:
                if flag in line:
                    N = int(line.split(',')[4])
        final_df = addRPKM(genome, sample, count_df, N)
        out_file = "{}_{}.csv".format(genome, sample)
        out_dir = config.get('output', 'path')
        out_path = os.path.join(out_dir, out_file)
        final_df.to_csv(out_path, index_label=genome)

def checkGETHOMS(genome, ref_genomes, GETHOMS):
    genome_list = os.path.join(GETHOMS, "input_order.txt")
    with open(genome_list, "r") as gl:
        input = gl.read()
    return all([g in input for g in ([genome] + ref_genomes)])



def isCrossRef(matrix):
    if os.path.isfile(matrix.split(".")[0] + "_crossRef.csv"):
        return matrix.split(".")[0] + "_crossRef.csv"
    return False

def cft073():

    gbk ="/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/" \
         "2018-01-16-sync-genomes/data/" \
         "CFT073/GCA_000007445.1_ASM744v1_genomic.gbff"
    genome = "CFT073"
    return gbk, genome

def MG1655():
    gbk = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/" \
                 "2018-01-16-sync-genomes/data/" \
                  "MG1655/GCA_000005845.2_ASM584v2_genomic.gbff"
    genome = "MG1655"
    return gbk, genome

def get_orthologs(pgm, genome, ref_genomes):
    # Given a genome and gene_id provide gene id of ortholog from a reference genome
    # Relies on get_homologues output
    #print("Index Genome: {}".format(index_genome))
    # Dropping all genes not found in reference, or found in multiple copies
    #print(pgm.shape)
    pgm = pgm[pgm[genome] != "PARALOGS"]
    pgm.dropna(subset=[genome], inplace=True)
    pgm.set_index(genome, inplace=True)
    result = pgm[ref_genomes]
    #result.to_csv("test_out/HMO1.csv")
    #print(result[:3])
    return result

# Actually want to pass in genome, condition and seq run(seq1, seq2 or nothing)
def addCounts(genome, sample, ortho_df, count_path, gff_path):
    # Read in counts
    co = pd.read_csv(count_path,sep ="\t", index_col=0, header=None)
    co.columns = ["{}_{}_counts".format(genome, sample)]
    todrop = [ind for ind in list(co.index) if "__" in ind]
    co.drop(todrop, inplace=True)
    out = co.join(ortho_df)
    gene_len = {}
    with open (gff_path, "r") as gff:
        for line in gff:
            if "CDS" not in line:
                continue
            try:
                gene = line.split(";")[0].split("gene_id=")[1]
                gene_len[gene] = abs(int(line.split()[4])-int(line.split()[3]))
            except:
                pass
    out["gene_length"] = pd.Series(gene_len)
    #print(out.head)
    return out

def addRPKM(genome, sample, count_df, N):

    count_df["{}_{}_RPKM".format(genome, sample)] = round((10**9*count_df["{}_{}_counts".format(genome, sample)])/(N*count_df["gene_length"]),2)
    return count_df


def combineCounts(samples, config):
    GETHOMS = config.get("GETHOMS", "path")
    path = os.path.join(GETHOMS, "intersect_pan{}".format(config.get("GETHOMS", "suffix")))
    assert os.path.isdir(path)
    matrix = os.path.join(path, "pangenome_matrix{}.tab".format(config.get("GETHOMS", "matrix_suffix")))
    crossRef = isCrossRef(matrix)
    pgm = pd.read_csv(crossRef)
    counts_df = pgm.copy()
    for counted in samples:
        genome = counted.split("_")[0]
        print(genome)
        path = os.path.join(config.get("output", "path"), "{}.csv".format(counted))
        counts = pd.read_csv(path)
        counts.drop(['CFT073', 'MG1655', 'gene_length'], axis = 1, inplace=True)
        print(counts[0:3])
        counts_df = counts_df.merge(counts, how='outer', on=genome)
        print(counts_df[0:3])
    out_dir = config.get("output", "path")
    out_file = os.path.join(out_dir, "{}_combined_counts.csv".format(datetime.datetime.now().strftime("%Y-%m-%d")))
    counts_df.to_csv(out_file)


def combineCustomCounts(samples, config):
    counts_df = pd.DataFrame()
    for counted in samples:
        genome = counted.split("_")[0]
        path = os.path.join(config.get("output", "path"), "{}.csv".format(counted))
        counts = pd.read_csv(path, index_col=0)
        counts.drop(['CFT073', 'MG1655', 'gene_length'], axis = 1, inplace=True)
        print(counts[0:3])
        if counts_df.empty:
            counts_df = counts.copy()
        else:
            counts_df =pd.merge(counts_df, counts, left_index=True, right_index=True)
        print(counts_df[0:3])
    out_dir = config.get("output", "path")
    out_file = os.path.join(out_dir, "{}_combined_custom_counts.csv".format(datetime.datetime.now().strftime("%Y-%m-%d")))
    counts_df.to_csv(out_file)




    """
    Calculate RPKM
 	RPKM = (10^9 * C)/(N * L), with

    # C = Number of reads mapped to a gene
    # N = Total mapped reads in the experiment
    # L = gene length in base-pairs for a gene

    """
