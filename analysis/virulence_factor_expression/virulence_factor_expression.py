import argparse
import datetime as dt
import os
import pandas as pd
import sys

sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/methods')

import keggAPI
import helpers
import run_blast

TODAY = dt.datetime.today().strftime("%Y-%m-%d")


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-vf", dest="virulence_file", required=False)
    parser.add_argument("-o", dest="outdir", required=False)
    parser.add_argument("-a", dest="analysis", help="mf or blast or both", required=False)
    parser.add_argument("-mf", dest="multifasta", required=False)
    return parser


def make_presence_absence_matrix(gene_to_prokka, config_dict, output_directory,
                                 genome_set="all"):

    strains = config_dict["genomes"][genome_set].split()
    pa_matrix = {}
    pa_matrix_prokka = {}
    for gene, prokka_list in gene_to_prokka.items():
        pa_matrix[gene] = {st: 0 for st in strains}
        pa_matrix_prokka[gene] = {st: '' for st in strains}
        for genome_prokka_tuple in prokka_list:
            prokka = genome_prokka_tuple[1]
            if "PROKKA" not in prokka:
                continue
            pa_matrix[gene][genome_prokka_tuple[0]] += 1
            # for paralogs this will only keep one of them
            pa_matrix_prokka[gene][genome_prokka_tuple[0]] = prokka
            # Will look at expression of representative
    matrix_file = os.path.join(output_directory,
                               "{}_presence_absence_matrix.csv".format(TODAY))
    matrix_prokka_file = os.path.join(output_directory,
                                      "{}_presence_absence_matrix_prokka.csv".format(TODAY))
    matrix = pd.DataFrame(pa_matrix).T
    matrix_prokka = pd.DataFrame(pa_matrix_prokka).T
    matrix.to_csv(matrix_file)
    matrix_prokka.to_csv(matrix_prokka_file)
    return matrix_file


def compare_virulence_gene_expression(virulence_factor_file="", output_directory="",
                                      config="", blast="nt", multi_fasta="",
                                      blast_run=True, clean_up=False, today=TODAY):
    if not config:
        config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")
    config_dict = helpers.process_config(config)
    if not output_directory:
        output_directory = config_dict["output_directory"]["path"]
    # 1 Get nucleotide sequence for each virulence factor
    if blast_run and not multi_fasta:
        if not virulence_factor_file:
            virulence_factor_file = config_dict["virulence_genes"]["path"]
        multi_fasta = keggAPI.get_genes(virulence_factor_file, '', output_directory, blast)

    # 2 For each gene run blast
    if blast_run:
        blast_bin = config_dict["blast"]["path"]
        genome_folder = config_dict["blast"]["genome_folder"]  # this is where clinical strains are
        db_prefix = config_dict["blast"]["db_prefix"]
        run_blast.run_nucleotide_blast(output_directory, blast_bin,
                                       db_prefix, genome_folder,
                                       multi_fasta, today=TODAY)

    # 3 Process Blast output
    gene_in_genomes, gene_identities, gene_coverage = run_blast.process_blast_output(output_directory)

    # 4 Identify Prokkas
    gff_folder = config_dict["gff_folder"]["path"]
    gene_to_prokka = run_blast.find_all_overlaps(gene_in_genomes, gff_folder)
    # I need this to be written out

    # 5 Build p/a matrix
    make_presence_absence_matrix(gene_to_prokka, config_dict,
                                 output_directory, genome_set="all")

    # 6 write out mean identity and coverage
    identity_file = os.path.join(output_directory, (today + "_identity_coverage.csv"))
    with open(identity_file, "w") as i_f:
        i_f.write('"",identity,coverage\n')
        for gene in gene_identities.keys():
            i_f.write('{},{},{}\n'.format(gene, gene_identities[gene], gene_coverage[gene]))

    # 7 Clean up: remove all the blast output files, etc


if __name__ == "__main__":

    args = get_args().parse_args()
    if args.analysis == "blast" or args.analysis == "both":
        blast_run = True
    else:
        blast_run = False
    compare_virulence_gene_expression(virulence_factor_file=args.virulence_file,
                                      output_directory=args.outdir,
                                      config="", blast="nt", multi_fasta=args.multifasta,
                                      blast_run=blast_run, clean_up=False)


