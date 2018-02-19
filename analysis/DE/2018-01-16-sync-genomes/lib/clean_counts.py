
"""


Given: a list of sample names (ex. HM01_UR, or HM07_UTI_seq1), create a csv file with counts, RPKMs and MG1655 and CFT073 orthologs
Also creates individual csv file for each sample

"""

import os
import configparser
import subprocess
import argparse
from modules import get_gene_id
from modules import generating_OGC_crossref as cr
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


def orthologFinder():





#config = configparser.ConfigParser()
#config.read("lib/human_config")
GETHOMS = config.get("GETHOMS", "path") # dir
HTSEQ = config.get("HTSEQ", "path") # dir
FLAG = config.get("FLAG", "path") # dir

if not os.path.isfile(os.path.join(FLAG, "flagstat_summary.txt")):
    pf.parseFlagstat(config)
assert os.path.isfile(os.path.join(FLAG, "flagstat_summary.txt"))


GFF = config.get("GFF", "path") # dir or file check
OUTDIR = config.get("output", "path") # dir


def clean_counts(samples, config, out_dir=OUTDIR): # This will only work as long as outdir is specified in the config file
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
    #samples = generateSampleNames(genomes = ["HM01", "HM06", 'HM07', 'HM57'], sample = ["UR", "UTI"], reseqed = False)
    samples = generateSampleNames()
    clean_counts(samples, config)

    #######################
    args = parser().parse_args()
    config = configparser.ConfigParser()
    if args.config:
        config.read(args.config)
    else:
        config.read(os.path.dirname(os.path.abspath(__file__)))







