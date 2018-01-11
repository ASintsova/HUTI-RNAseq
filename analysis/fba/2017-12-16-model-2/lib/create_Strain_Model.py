from __future__ import print_function
import cobra
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
import os
#import argparse
import re
import pandas

BASE_MODEL_PATH="/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data/ref/mg1655_2011_model.xml"
BASE_MODEL_NAME=os.path.basename(BASE_MODEL_PATH).split(".")[0]
FOLDER = '/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data'

#def parser():
#    parser = argparse.ArgumentParser()
#    parser.add_argument("-f", "--ogc_file",  help= "Orthologous Groups Clusters File", required = True)
#    return parser

def keepEssentials(model, gene_list):
    growth_rate_factor = model.optimize().objective_value/10.0 # This is arbitrary
    #for each gene in gene list check if knocking out breaks the model
    deletion_results = single_gene_deletion(model, gene_list)
    to_keep = deletion_results.flux > growth_rate_factor
    new_genes_to_remove = list(deletion_results[to_keep].index.values)
    return new_genes_to_remove

def generateStrainModel(ogc_file):
    # Get all the genes found in clinical UTI
    strain = os.path.basename(ogc_file).split("OGC")[0]
    print("Running generateStrainModel for {}".format(strain.strip("_")))
    print("Getting all the genes from {}".format(os.path.basename(ogc_file)))
    bnum = re.compile(r'b\d\d\d\d')
    with open(ogc_file, "r") as fh:
        text = fh.read()
        strain_genes = bnum.findall(text)
    print("{} genes found. First 10 are: {}".format(len(strain_genes), ", ".join(strain_genes[:10])))

    # Load base model
    print("Loading base model: {}".format(BASE_MODEL_NAME))
    base = cobra.io.read_sbml_model(BASE_MODEL_PATH)
    new_model_name= strain + base.id
    model_genes = list(base.genes)
    original_num_reactions = len(base.reactions)
    original_num_genes = len(base.genes)

    # Find genes not found in clinical strain
    print("Looking for genes not found in the clinical strain")
    genes_to_remove = [g for g in model_genes if str(g) not in strain_genes]
    print("Found {} genes to remove".format(len(genes_to_remove)))

    # Check if any of them are essential
    print("Checking if any of them are 'essential'")
    genes_to_remove_checked_for_phenotypes = keepEssentials(base, genes_to_remove)
    print("Found {} 'essential' genes".format(len(genes_to_remove) - len(genes_to_remove_checked_for_phenotypes)))
    # Make a list of essential genes not found in the clinical strain
    essentials = [str(g) for g in genes_to_remove if str(g) not in genes_to_remove_checked_for_phenotypes]
    #Can then see if they are actually missing or if it's a faulty annotation
    print("Essential genes: ", str(essentials))

    # Delete appropriate reactions from base model
    print ("Deleting genes and reactions")
    if len(genes_to_remove_checked_for_phenotypes)>0:
        cobra.manipulation.delete.remove_genes(base, genes_to_remove_checked_for_phenotypes, remove_reactions=True)

    new_num_reactions = len(base.reactions)
    new_num_genes = len(base.genes)
    print ("There were {} genes to remove, actually removed {}".format(len(genes_to_remove_checked_for_phenotypes), (original_num_genes - new_num_genes)))
    print("There were {} reactions, now there are {}".format(original_num_reactions, new_num_reactions))

    # Write out the new model
    print("Writing out new model")
    base.id = new_model_name
    output_dir = os.path.join(os.path.dirname(ogc_file), (strain +"model"))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    model_file = os.path.join(output_dir, (new_model_name+".xml"))
    cobra.io.write_sbml_model(base, model_file)

    essential_genes_out_file = os.path.join(output_dir, (strain + "retained_NOT_FOUND.txt"))
    deleted_genes_out_file = os.path.join(output_dir, (strain + "deleted_NOT_FOUND.txt"))
    with open(essential_genes_out_file, "w") as eg:
        eg.write("\n".join(essentials))
    with open (deleted_genes_out_file, "w") as dg:
        dg.write("\n".join(genes_to_remove_checked_for_phenotypes))

#args = parser().parse_args()
#ogc_file = args.ogc_file

def findOGCFiles(folder):
    ogc_file_list =[]
    for root, dirs,files in os.walk(folder, topdown =True):
        for fname in files:
            if "OGC.tab" in fname:
                ogc_file_list.append(os.path.join(root, fname))
    return ogc_file_list

for ogc_file in findOGCFiles(FOLDER):
    generateStrainModel(ogc_file)
