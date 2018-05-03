import os
import pandas as pd
import datetime as dt
from modules import keggAPI

def deletedGenes(data_folder = "./data", output_dir="./data"):

    today = dt.datetime.today().strftime("%Y-%m-%d")
    deleted_genes = {}
    deleted_reactions = {}
    essential_genes_NOT_FOUND = {}

    for root, dirs,files in os.walk(data_folder, topdown =True):
        for fh in files:
            if "deleted_NOT_FOUND" in fh:
                with open(os.path.join(root, fh)) as fi:
                    for line in fi:
                        line = line.rstrip("\n")
                        deleted_genes[line] = deleted_genes.get(line, 0) + 1
            elif "retained_NOT_FOUND" in fh:
                with open(os.path.join(root, fh)) as fi:
                    for line in fi:
                        line = line.rstrip("\n")
                        essential_genes_NOT_FOUND[line] = essential_genes_NOT_FOUND.get(line, 0) + 1
            elif "deleted_REACTIONS" in fh:
                with open(os.path.join(root, fh)) as fi:
                    for line in fi:
                        line = line.rstrip("\n")
                        deleted_reactions[line] = deleted_reactions.get(line,0) + 1


    bnums = pd.Series(deleted_genes,index=deleted_genes.keys())
    bnums = bnums.sort_values(ascending=False)
    bnums.to_csv(os.path.join(output_dir, "{}-deleted_genes.csv".format(today)))

    essent_bnums = pd.Series(essential_genes_NOT_FOUND, index=essential_genes_NOT_FOUND.keys())
    essent_bnums = essent_bnums.sort_values(ascending=False)
    essent_bnums.to_csv(os.path.join(output_dir, "{}-retained_NOT_FOUND.csv".format(today)))
    reactions = pd.Series(deleted_reactions, index=deleted_reactions.keys())
    reactions = reactions.sort_values(ascending=False)
    reactions.to_csv(os.path.join(output_dir, "{}-deleted_reactions.csv".format(today)))
    return (os.path.join(output_dir, "{}-deleted_genes.csv".format(today)),
        os.path.join(output_dir, "{}-retained_NOT_FOUND.csv".format(today)),
        os.path.join(output_dir, "{}-deleted_reactions.csv".format(today)))


def getInfoOnDeletedGenes(deleted_genes_csv, strain_cutoff, genome="eco"):
    today = dt.datetime.today().strftime("%Y-%m-%d")
    df = pd.read_csv(deleted_genes_csv, index_col=0, names=["strains"])
    df2 = df[df.strains == strain_cutoff]
    df2 = df2.sort_index()
    gene_call_list = ["{}:{}".format(genome, bnum) for bnum in df2.index]
    print(gene_call_list)
    return gene_call_list, "data/{}-gene-info-deleted-from-{}-strains".format(today, strain_cutoff)


x, _, _ = deletedGenes()

gene_call_list, output_prefix = getInfoOnDeletedGenes(x, 13)
keggAPI.getgeneInfo(gene_call_list, output_prefix)