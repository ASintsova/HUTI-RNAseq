"""
Given two gbk files and a base model (corresponding to one of the gbk files
create a strain specific model

"""

from __future__ import print_function
import cobra
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
import cobra.test
import os
import re
import shlex
import subprocess
from modules import mapOGCs
import glob
import argparse
import configparser

##############################################################################################
"""

Set up: data folder which contains a subfolder for each model to be build containing that strains gbk + model gbk

"""
##############################################################################################
# Get list of strains from user

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--strains', nargs='+', help = 'list of strains to generate models',
                            type=str, required=True)
    return parser


##############################################################################################
# Run BDBH on the two genebank files

def runGH(data_dir, strain,  GH_path, S, C):
    cwd = os.getcwd()
    os.chdir(data_dir)
    if os.path.isdir(strain):
        gbk_dir = strain
        cmd_str = "perl {} -d {} -S {} -C {}".format(os.path.join(GH_path, "get_homologues.pl"), gbk_dir, S, C)
        print(cmd_str)
        cmd = shlex.split(cmd_str)
        subprocess.Popen(cmd)
    else:
        print("No strain directory found")
    os.chdir(cwd)

##############################################################################################

# Generate ortholog map
def getClusterMap(bnum_file, hom_dir):

    fasta_dir = glob.glob(hom_dir+"/*BDBH*")[0]
    cL = glob.glob(hom_dir+"/*.cluster_list")[0]
    out = os.path.join(hom_dir, os.path.basename(hom_dir).split("_")[0]+"_OGC.tab")
    mapOGCs.mapOGCs(cL, fasta_dir, bnum_file, out)

# Again will probably want to improve on hardcoding directories in

#######################################################################################

def keepEssentials(model, gene_list):
    # Want this to be 5% of normal
    growth_rate_factor = model.optimize().objective_value*0.05

    # For each gene in gene list check if knocking out breaks the model

    deletion_results = single_gene_deletion(model, gene_list)
    to_keep = deletion_results.growth > growth_rate_factor
    new_genes_to_remove = [ge for gene in list(deletion_results[to_keep].index) for ge in gene]
    return new_genes_to_remove

def essentialReactions(essentials, base_model_path):

    """
    :param essentials: list of essential genes
    :param base_model_path: path to the base model
    :return: set of essential reactions that could not be deleted by has not genes info to support them

    """
    if base_model_path.endswith("ml"):
        base = cobra.io.read_sbml_model(base_model_path)
    else:
        base = cobra.io.load_json_model(base_model_path)
    orx = [r.id for r in base.reactions]
    cobra.manipulation.delete.remove_genes(base, essentials, remove_reactions=True)
    nrx = [r.id for r in base.reactions]
    return set(orx) - set(nrx)


def generateStrainModel(ogc_file, base_model_path, format):
    """
    :param ogc_file: tab file mapping orthologs between strain of interest and base model genes
    :param base_model_path: path to the base model
    :return: writes out new strain specific model, genes that were deleted, genes that could not be deleted, and reactions that were deleted

    """
    base_model_name = os.path.basename(base_model_path).split(".")[0]
    # Get all the genes found in clinical UTI
    strain = os.path.basename(ogc_file).split("_OGC")[0]
    print("Running generateStrainModel for {}".format(strain))
    print("Getting all the genes from {}".format(os.path.basename(ogc_file)))
    bnum = re.compile(r'b\d\d\d\d')
    with open(ogc_file, "r") as fh:
        text = fh.read()
        strain_genes = bnum.findall(text) + ['s0001']
    print("{} genes found.".format(len(strain_genes)))

    # Load base model
    print("Loading base model: {}".format(base_model_name))
    if format == 'JSON':
        base = cobra.io.load_json_model(base_model_path)
    else:
        try:
            base = cobra.io.read_sbml_model(base_model_path)
        except:
            print("Could not determine model format")

    original_growth = base.optimize().objective_value
    original_reactions = [r.id for r in base.reactions]

    new_model_name= strain +"_" + base.id
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
    if len(genes_to_remove_checked_for_phenotypes) > 0:
        cobra.manipulation.delete.remove_genes(base, genes_to_remove_checked_for_phenotypes, remove_reactions=True)
    new_num_reactions = len(base.reactions)
    new_num_genes = len(base.genes)
    print ("There were {} genes to remove, actually removed {}".format(len(genes_to_remove_checked_for_phenotypes), (original_num_genes - new_num_genes)))
    print("There were {} reactions, now there are {}".format(original_num_reactions, new_num_reactions))
    new_growth = base.optimize().objective_value
    new_reactions = [r.id for r in base.reactions]

    deleted_reactions = set(original_reactions) - set(new_reactions)
    print("Deleted reactions: {}".format(deleted_reactions))
    print("Used to grow at rate of {}. Now grows at rate of {}".format(original_growth, new_growth))
    # Write out the new model
    print("Writing out new model")
    base.id = new_model_name
    output_dir = os.path.join(os.path.dirname(ogc_file), (strain +"_model"))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    model_file = os.path.join(output_dir, (new_model_name+".xml"))
    cobra.io.write_sbml_model(base, model_file)

    essential_genes_out_file = os.path.join(output_dir, (strain + "_retained_NOT_FOUND.txt"))
    deleted_genes_out_file = os.path.join(output_dir, (strain + "_deleted_NOT_FOUND.txt"))
    deleted_reactions_out_file = os.path.join(output_dir, (strain + "_deleted_REACTIONS.txt"))

    with open(essential_genes_out_file, "w") as eg:
        eg.write("\n".join(essentials))
    with open (deleted_genes_out_file, "w") as dg:
        dg.write("\n".join(genes_to_remove_checked_for_phenotypes))
    with open (deleted_reactions_out_file, "w") as dr:
        dr.write("\n".join(deleted_reactions))
    print(essentialReactions(essentials, base_model_path))


def pipeline():
    config = configparser.ConfigParser()
    config.read("lib/config")
    data_dir = config.get("paths", "data")

    args = parser().parse_args()

    # Get list of strains
    strains = args.strains
    for st in strains:
        print("Strain: {}".format(st))
        # Check for get_homologues output
        if not os.path.isdir(os.path.join(data_dir, "{}_homologues".format(st))):
            gh_path = config.get("get_homologues", "path")
            S = config.get("get_homologues", "S")
            C = config.get("get_homologues", "C")
            runGH(data_dir, st, gh_path, S, C)
        else:
            print("Found blast output")
        # Check for OGC file

        if not os.path.isfile(os.path.join(data_dir,
                             "{}_homologues".format(st),
                            "{}_OGC.tab".format(st))):
            print("No OGC file")

            # Generate OGC file:
            bnum_file = config.get("paths", "bnum_file")
            hom_dir = os.path.join(data_dir,"{}_homologues".format(st))
            getClusterMap(bnum_file, hom_dir)
        else:
            print("Found OGC file")
        ogc_file = os.path.join(data_dir,
                             "{}_homologues".format(st),
                            "{}_OGC.tab".format(st))

        base_model_path = config.get("base_model", "path")
        fmt = config.get("base_model", "format")
        generateStrainModel(ogc_file, base_model_path, fmt)

if __name__ == "__main__":

    pipeline()
