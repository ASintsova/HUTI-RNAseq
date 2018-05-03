
"""
Appropriate file naming formats are outlined in the config file as well

Given: a list of sample names (ex. HM01_UR, or HM07_UTI_seq1), create a csv file with counts, RPKMs and MG1655 and CFT073 orthologs
Also creates individual csv file for each sample


1. Run get_homologues on genomes of interest, specific configuration to be found in the config file
2. Using output from get_homolouges generate a crossReference (module crossRefOGC) matching orthologs between different genome annotations

3.

"""

import os
import configparser
import subprocess
import argparse
from modules import get_gene_id
from modules import getOGCCounts as gc
from modules import parsing_gbk as pg
from itertools import product
from modules import parse_flagstats as pf



######################################################################################
#### Helper function to generate all possible, or any combination of sample names for my clinical strains.
#### Should probably be moved somewhere else

ALL = ["HM01", "HM03","HM06", "HM07",
              "HM14", "HM17", "HM43", "HM54", "HM56",
              "HM57", "HM60", 'HM68',"HM66", "HM86"]
COND = ["UR", "UTI"]

def generateSampleNames(genomes = ALL, sample = COND, reseqed = True):
    reseq = ["HM06", "HM07", "HM43", "HM57", "HM60", 'HM68']

    sample_names = list(product(genomes, sample))
    if reseqed:
        attr = ["UTI_seq1", "UTI_seq2"]
        additional_samples = [ge for ge in genomes if ge in reseq]
        additional_names = list(product(additional_samples, attr))
        sample_names += additional_names
    names = ["_".join(s) for s in sample_names]
    return names

#######################################################################################

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=False)
    parser.add_argument("-o", "--out_dir", required=True)
    parser.add_argument("-s", "--samples", nargs='+', help="List of samples", required=False)



#def orthologFinder():
    #this is where running of get_homologues would go


########################################################


########################################################

# 1. Check if CrossRef has been generated

def isCrossRef(matrix):
    if os.path.isfile(matrix.split(".")[0] + "_crossRef.csv"):
        return matrix.split(".")[0] + "_crossRef.csv"
    return False

# 2. Check that gethomologues output contains all the genomes you're interested in

def checkGETHOMS(genome, ref_genomes, GETHOMS):
    genome_list = os.path.join(GETHOMS, "input_order.txt")
    with open(genome_list, "r") as gl:
        input = gl.read()
    return all([g in input for g in ([genome] + ref_genomes)])


def getOGCCounts(samples, config, out_dir):

    subprocess.call(["mkdir", "-p", out_dir])
    assert os.path.isdir(out_dir)

    for sample in samples:
        genome = sample.split("_",1)[0]
        #print(genome)
        cond = sample.split("_",1)[1]
        #print(cond)
        get_gene_id.get_gene_id(genome, cond, ['CFT073', 'MG1655'], config)

    #get_gene_id.combineCounts(samples, config)
    get_gene_id.combineCustomCounts(samples, config)
    ### DEBUGGING: can't open crossRef file


if __name__ == "__main__":
    # Get info from config parser


    args = parser().parse_args()
    config = configparser.ConfigParser()
    if args.config:
        config.read(args.config)
    else:
        config.read(os.path.dirname(os.path.abspath(__file__)))

    GETHOMS = config.get("GETHOMS", "path")  # dir
    HTSEQ = config.get("HTSEQ", "path")  # dir
    FLAG = config.get("FLAG", "path")  # flagstat_summary.txt
    GFF = config.get("GFF", "path")  # dir or file check


    if args.samples:
        samples = [sam for sam in args.input]
    else:
        samples = samples = generateSampleNames()










