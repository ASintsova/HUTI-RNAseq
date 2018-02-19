# Needs file path to pangenome matrix, and to the directory of 'intersect_pan_CM_t0' with *faa files
# Will return a csv file with all gene_ids cross-referenced
# This will have protein_id for MG1655 and CFT073, need to parse gbk to parse out locus tags
# Missing genes will be represented by 0, multiple copies by appropriate integer


import pandas as pd
import os
import re
import numpy as np

test = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes/test_out/test_matrix.tab"
path = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes" \
       "/data/genomes_homologues/intersect_pan_CM_t0"

def crossRefFromPanGenome(matrix, pangene_seq_dir):
    # Read in pangenome matrix
    pgm = pd.read_csv(matrix, sep="\t", index_col=0).T
    # Remove 'gbf' from column names
    col_names = {c: c.split(".")[0] for c in list(pgm.columns)}
    pgm.rename(columns=col_names, inplace=True)
    pgm.columns.name = 'gh_id'
    # Go through each cell, and replace 1 with appropriate gene_ids
    for index, row in pgm.iterrows():
        for col in list(pgm.columns):
            if row[col] == 1:
                if "HM0" in col:
                    genome_re = re.compile(col+"_|" + col.replace("0", "")+"_")
                else:
                    genome_re = re.compile(col)

                with open(os.path.join(pangene_seq_dir, index), "r") as fh:
                    for line in fh:
                        if genome_re.search(line):
                            # get ">ID:" put it instead of one.
                            gene_id = line.split("|")[0].lstrip(">ID:").rstrip()
                            pgm.loc[index, col] = gene_id
            elif row[col] > 1:
                pgm.loc[index, col] = "PARALOGS"
            elif row[col] == 0:
                pgm.loc[index, col] = np.nan
            else:
                continue
                # Or replace with NaN?
    out_file = os.path.join(pangene_seq_dir, (matrix.split(".")[0]+"_crossRef.csv"))
    pgm.to_csv(out_file)
    return out_file

#crossRefFromPanGenome(test, path)