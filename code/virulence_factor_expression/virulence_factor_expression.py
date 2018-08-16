import datetime as dt
import glob
import os
import pandas as pd
import sys

sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/code/methods')

import keggAPI
import helpers
import run_blast


TODAY = dt.datetime.today().strftime("%Y-%m-%d")


def get_rpkms_for_prokka(rpkm_folder, genome,  prokka, condition):

    filename = glob.glob(os.path.join(rpkm_folder, "*{}_{}*".format(genome, condition)))
    if not filename:
        genome = "HM0"+ genome.split("HM")[1]
        filename = glob.glob(os.path.join(rpkm_folder, "*{}_{}*".format(genome, condition)))
    with open(filename[0], "r") as fh:
        for line in fh:
            if prokka in line:
                rpkm = line.split(",")[1].strip()
                return float(rpkm)


def make_presence_absence_matrix(gene_to_prokka, config_dict, output_directory,
                                 genome_set="all"):
    Genome = config_dict["virulence_genes"]["genome"]
    strains = config_dict["genomes"][genome_set].split()
    pa_matrix = {}
    pa_matrix_prokka = {}
    for gene, prokka_list in gene_to_prokka.items():
        pa_matrix[gene] = {st:0 for st in strains}
        pa_matrix_prokka[gene] = {st:0 for st in strains}
        for st in prokka_list:
            prokka = st[1]
            if "PROKKA" not in prokka:
                continue
            pa_matrix[gene][st[0]] = 1
            pa_matrix_prokka[gene][st[0]] = prokka
    matrix_file = os.path.join(output_directory,
                               "{}_{}_presence_absence_matrix.csv".format(TODAY, Genome))
    matrix_prokka_file = os.path.join(output_directory,
                               "{}_{}_presence_absence_matrix_prokka.csv".format(TODAY, Genome))
    matrix = pd.DataFrame(pa_matrix).T
    matrix_prokka = pd.DataFrame(pa_matrix_prokka).T
    matrix.to_csv(matrix_file)
    matrix_prokka.to_csv(matrix_prokka_file)
    return matrix_file


def get_rpkms_for_virulence_factor(gene_to_prokka, config_dict, output_directory='.', genome_set="all",
                                 conditions = ["UR", "UTI"]):

    rpkm_folder = config_dict["counts"]["rpkm_path"]
    rpkm_matrix = {}
    for gene, prokka_list in gene_to_prokka.items():
        rpkm_matrix[gene] = {}
        for cond in conditions:
            for st in prokka_list:
                genome = st[0]
                cond_name = "{}_{}".format(genome, cond)
                prokka = st[1]
                if "PROKKA" not in prokka:
                    continue
                rpkm = get_rpkms_for_prokka(rpkm_folder,
                                         genome, prokka,
                                         cond)
                rpkm_matrix[gene][cond_name] = rpkm

    matrix_file = os.path.join(output_directory,
                               "{}_RPKM_matrix.csv".format(TODAY))
    matrix = pd.DataFrame(rpkm_matrix).T
    matrix.to_csv(matrix_file)


def compare_virulence_gene_expression(config_dict, blast="nt"):

    # 1 Get nucleotide sequence for each virulence factor
    # todo refactor: make one virulence factor list, each line is genome, locus_tag
    gene_list_file = config_dict["virulence_genes"]["path"]
    #ref_genome = config_dict["virulence_genes"]["genome"]
    output_directory = config_dict["output_directory"]["path"]
    multi_fasta = keggAPI.get_genes(gene_list_file, output_directory, blast)  # todo make sure this works

    # 2 For each gene run blast
    blast_bin = config_dict["blast"]["path"]
    genome_folder = config_dict["blast"]["genome_folder"]  # this is where clinical strains are
    db_prefix = config_dict["blast"]["db_prefix"]

    run_blast.run_nucleotide_blast(output_directory, blast_bin,
                                   db_prefix, genome_folder,
                                   multi_fasta, today=TODAY)  # todo does this check/build db_index?

    # 3 Process Blast output
    gene_in_genomes = run_blast.process_blast_output(output_directory)

    # 4 Identify Prokkas
    gff_folder = config_dict["gff_folder"]["path"]
    gene_to_prokka = run_blast.find_all_overlaps(gene_in_genomes, gff_folder)
    # I need this to be written out

    # 5 Build p/a matrix

    pa_matrix_file = make_presence_absence_matrix(gene_to_prokka, config_dict, output_directory,
                                 genome_set="all")
    # 6 Build rpkm matrix

    get_rpkms_for_virulence_factor(gene_to_prokka, config_dict, output_directory, genome_set="all",
                                    conditions=["UR", "UTI"])

    # 7 Clean up: remove all the blast output files, etc





if __name__ == "__main__":
    config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")
    config_dict = helpers.process_config(config_file)
    compare_virulence_gene_expression(config_dict)


