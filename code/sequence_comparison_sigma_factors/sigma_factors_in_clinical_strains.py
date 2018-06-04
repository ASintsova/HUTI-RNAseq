import datetime as dt
import os
import sys
sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/code/methods')

import keggAPI
import helpers
import run_blast

TODAY = dt.datetime.today().strftime("%Y-%m-%d")

def get_CFT073_sigmas(blast = 'nt'):

    gene_list_file = config_dict["sigmas"]["path"]
    ref_genome = config_dict["sigmas"]["genome"]
    output_directory = config_dict["output_directory"]["path"]

    multi_fasta = keggAPI.get_genes(gene_list_file, ref_genome,
                        output_directory, blast )
    print(multi_fasta)


def find_clinical_strains_sigmas():
    blast_bin = config_dict["blast"]["path"]
    genome_folder = config_dict["blast"]["genome_folder"]
    db_prefix = config_dict["blast"]["db_prefix"]
    output_directory = config_dict["output_directory"]["path"]
    reference_sigma_fasta = config_dict['blast']['cft073_sigmas_nt']
    run_blast.run_nucleotide_blast(output_directory, blast_bin, db_prefix, genome_folder,
                                   reference_sigma_fasta, today=TODAY)

def process_blast_results():
    output_directory = config_dict["output_directory"]["path"]
    gene_in_genomes = run_blast.process_blast_output(output_directory)
    gff_folder = config_dict["gff_folder"]["path"]
    gene_to_prokka = run_blast.find_all_overlaps(gene_in_genomes, gff_folder)
    return gene_to_prokka



if __name__ == "__main__":
    config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")
    config_dict = helpers.process_config(config_file)
    process_blast_results()

