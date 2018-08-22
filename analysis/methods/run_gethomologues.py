"""

run get_homologues given a set of parameters

parameters I want to be able to change:

-C 50/60/70/75(default)
-S  1(default)/70/90
-e exclude paralogs

always to include
 -A
 -c


In the end want to check for # of singletons: go through all files in intersect_pan_...
and count how many files have only one entry, do a counting dictionary can then see distribution


Do I want to do this on smaller number of genomes? Yes, let's do 5: CFT073, HM01, HM07, HM57, HM66

Do I want to adapt this to cluster scenario? Not yet.
      for compare clusters also want to try -s (synteny) option and maybe -I
      (produce clusters with single-copy seqs from ALL taxa in file)
      -m is for pangenomic matrices and -T for the tree
        which I actually don't really need that

"""

import os
import subprocess
import shlex
import glob
import shutil


# get_homologs output directory will be the same as the one from which this script is run


def run_get_homologs(gh_bin, gbk_dir, I, C, S, exclude_paralogs, run_id):

    # Make a temp directory and go there
    if not os.path.exists(run_id):
        os.makedirs(run_id)
    os.chdir(run_id)

    # Copy gbk files specified in I into a new genome directory
    os.makedirs('genomes')
    with open(I, "r") as gf:
        for line in gf:
            print(line)
            file = os.path.join(gbk_dir, line.strip())
            shutil.copy(file, 'genomes')

    # Run gh with appropriate params
    # Now there's not I, and gbk dir is hard coded
    e = " -e" if exclude_paralogs else ''
    e_name = "1" if exclude_paralogs else "0"
    # Running COGtriangle algorithm
    gh_script = os.path.join(gh_bin, "get_homologues.pl")
    cmd1_str = "{} -d genomes -G -t 0 -C {} -S {} -A -c{}".format(gh_script, C, S, e)
    print(cmd1_str)
    cmd1 = shlex.split(cmd1_str)
    subprocess.call(cmd1)
    # Running OMCL algorithm
    cmd2_str = "{} -d genomes -M -t 0 -C {} -S {} -A -c{}".format(gh_script, C, S, e)
    cmd2 = shlex.split(cmd2_str)
    subprocess.call(cmd2)
    gh_out_dir = "genomes_homologues"

    cog_out = glob.glob(gh_out_dir+"/*algCOG_e{}_C{}_S{}_".format(e_name, str(C), str(S)))[0]
    omcl_out = glob.glob(gh_out_dir+"/*algOMCL_e{}_C{}_S{}_".format(e_name, str(C), str(S)))[0]

    if cog_out.split("alg")[0] == omcl_out.split("alg")[0]:
        print("Go!")
        # s = " -s" if synteny else ''

        # Run pangenome analysis
        cmd3_str = "{} -o {}/run_{}_pan_C{}_S{} " \
                   "-d {},{},  -t 0 -m -T".format(os.path.join(gh_bin, "compare_clusters.pl"),
                                                  gh_out_dir, run_id, str(C), str(S),
                                                  cog_out, omcl_out)
        print(cmd3_str)
        cmd3 = shlex.split(cmd3_str)
        print(cmd3)
        subprocess.call(cmd3)

    # Delete everything except for .tab and _list files
    for root, dirs, files in os.walk("genomes_homologues", topdown=False):
        for name in files:
            if '.tab' in name or '_list' in name:
                print(os.path.join(root, name))
            else:
                os.remove(os.path.join(root, name))
        for name in dirs:
            fold = os.path.join(root, name)
            print(fold)
            for f in os.listdir(fold):
                shutil.move(os.path.join(fold, f), root)
        for name in dirs:
            if not os.listdir(os.path.join(root, name)):
                os.removedirs(os.path.join(root, name))

import pandas as pd
import os
import re
import numpy as np


def crossRefFromPangenome(pangenome_matrix, OGC_dir): # only need one input really, replace with glob.glob

    """
    This function takes in the pangenome matrix and directory of OGC generated by get_homologues
    and returns a table of all orthologous genes. Ignores paralogs.

    Note to self: for MG1655 and CFT073 get_homologues gives protein_id not locus_tags, need to
    parse out locus_tags from gbk files.


    Missing genes will be represented by NaN, multiple copies by "PARALOGS"

    Generates csv file in the same directory as pangenome matrix ending with 'crossRef.csv'


    """

    # Read in the pangenome matrix
    pgm = pd.read_csv(pangenome_matrix, sep="\t", index_col=0).T

    # Remove 'gbf' from column names
    pgm.columns = [c.split(".")[0] for c in list(pgm.columns)]  # Columns are names of the genomes
    pgm.columns.name = ''

    # Go through each cell of pangenome matrix, and replace 1 with appropriate gene_ids
    for index, row in pgm.iterrows():  # Index is the name of orthologous cluster

        for col in list(pgm.columns):
            if row[col] == 1:
                if "HM0" in col:  # Some files are named HM07, some HM7, want to be able to find both
                    genome_re = re.compile(col+"_|" + col.replace("0", "")+"_")
                else:
                    genome_re = re.compile(col)

                with open(os.path.join(OGC_dir, index), "r") as fh:
                    for line in fh:
                        if genome_re.search(line):
                            # Get ">ID:" put it instead of one.
                            gene_id = line.split("|")[0].lstrip(">ID:").rstrip()
                            pgm.loc[index, col] = gene_id

            elif row[col] > 1:
                pgm.loc[index, col] = "PARALOGS"
            elif row[col] == 0:
                pgm.loc[index, col] = np.nan
            else:
                continue  # Hopefully this doesn't happen

    out_file = os.path.join(OGC_dir, (pangenome_matrix.split(".")[0]+"_crossRef.csv"))
    pgm.to_csv(out_file) # maybe add the date
    return out_file  # Returns the path to the file

if __name__ == "__main__":
    pangenome_matrix = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/get_homologs_output/C50_S90_e0_/run_C50_S90_e0__pan_C50_S90/pangenome_matrix_t0.tab"
    OGC_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/get_homologs_output/C50_S90_e0_/run_C50_S90_e0__pan_C50_S90"
    crossRefFromPangenome(pangenome_matrix, OGC_dir)