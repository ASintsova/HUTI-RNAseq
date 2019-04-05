import gseapy
import networkx as nx
import sys
import pandas as pd




import utilities as ut
import interacting_with_regulondb as ir

config_dict = ut.process_config("config")


# Strain Info

strain_info = pd.read_csv(strain_info_file, index_col=0)

# Regulon info
regulon_csv = config_dict["out_dir"]["regulon_csv"] # This was created using interacting_with_regulondb.py
regulon_gmt = config_dict["out_dir"]["regulon_gmt"]

# Counts/TPMs
deseq2_cnts_file = config_dict["data"]["deseq2_counts"]
sign_de_file = config_dict["data"]["sign_de_file"]



regulon = pd.read_csv(regulon_csv, index_col=0)
regulon.head()

# Create gene set file
###
if os.path.isfile(regulon_gmt):
    os.remove(regulon_gmt)
for group_name, group  in regulon.groupby(["regulator_function", "regulator"]):
    with open(regulon_gmt, "a") as fo:
        fo.write("{}_{}\tNA\t{}\n".format(group_name[1],group_name[0],"\t".join(group.regulated_bnum.values)))


# Using DESeq2 normalized counts
ds_cnts =  pd.read_csv(deseq2_cnts_file, index_col=0)
print(ds_cnts.iloc[:, 0:4].head())
print(ds_cnts.shape)

# Setting up conditions vector
cls = ["UR","UTI"]*14


# Function to run gsea and edit resulting df

def run_gsea(data, gene_set, cls, outdir, fdr=0.05, file_name="test_run.csv"):
    gsea_run = gseapy.gsea(data=data, gene_sets=gene_set, cls=cls, outdir=outdir,
                           min_size=10, method="signal_to_noise")
    result = pd.DataFrame(gsea_run.results).T
    result = result[["es", "nes", "pval", "fdr", "geneset_size", "matched_size",
                     "genes"]].sort_values('fdr')
    results_short = (result[["nes", "fdr", "geneset_size", "matched_size", "genes"]][
                         result.fdr < fdr].sort_values(["nes", "fdr"]))
    results_short.columns = ["NES", "FDR", "Regulon Size", "Matched Size", "Genes"]
    results_short.to_csv(os.path.join(outdir, file_name))
    return results_short


def gene_name_to_protein_name(name):
    return name[:3].capitalize() + name[3:].capitalize()


###
gene_set = regulon_gmt
ga = run_gsea(ds_cnts, gene_set, cls, results_dir, file_name='final_gsea_analysis.csv')
