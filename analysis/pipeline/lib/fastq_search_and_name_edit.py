import os
import pandas as pd
import subprocess

FASTQ_DIR = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/fastq"
READS_DIR = "/Volumes/Passport_I/ANNA/clinical_strains_rnaseq/data/reads"
CSV_FILE = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/core_sample_numbers.csv"

sample_nums = pd.read_csv(CSV_FILE, index_col=0)

for root, dirs,files in os.walk(READS_DIR, topdown =True):

    for f_name in files:
        if ".fastq" in f_name:
            file_path = os.path.abspath(os.path.join(root, f_name))

            core_tag = f_name.split("_")[0]
            #print(core_tag)
            try:
                new_name_tag = sample_nums.loc[int(core_tag)]["Sample #"]
                new_name = new_name_tag +"_" +"_".join(f_name.split("_")[1:])
                new_path = os.path.join(FASTQ_DIR, new_name)
                subprocess.call(["cp", file_path, new_path])
            except:
                pass

