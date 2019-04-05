"""

This script requires config file to contain the following:

- valid path to locally installed BLAST
- valid path to folder with fasta files for all genomes of interest
- valid path to folder with gff files for all genomes of interest

TODO:

- Add documentation
- Add fasta folder as an argument
- Add gff folder as an argument

- Add % identity/coverage as arguments
- Write a simple example
- Refactor
- Add clean up

"""



import argparse
import datetime as dt
import numpy as np
import os
import pandas as pd
import pybedtools
import subprocess
import shlex
import sys


from kegg_gene import Gene, GeneSet
import utilities as ut


TODAY = dt.datetime.today().strftime("%Y-%m-%d")


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-vf", dest="virulence_file", required=False)
    parser.add_argument("-o", dest="outdir", required=False)
    parser.add_argument("-a", dest="analysis", help="mf or blast or both", required=False)
    parser.add_argument("-mf", dest="multifasta", required=False)
    return parser

def make_nucleotide_blast_db(genome_folder, output_directory, db_prefix,
                             blast_bin):

    """
    Checks if db exists, if not makes one.

    :param genome_folder: folder conatining fasta genomes to be used for custom db
    :param output_directory: directory where to put db.
    :param db_prefix: db name
    :param blast_bin: location of blast tool on your computer
    :return: db_path - path to the file containing
    concatenation of all genomes in the db

    """
    genomes = [os.path.join(genome_folder, g) for g in os.listdir(genome_folder)]
    db_path = os.path.join(output_directory, (db_prefix + ".fsa"))
    if not os.path.isfile(db_path):
        with open(db_path, 'w') as outfile:
            for genome in genomes:
                with open(genome) as infile:
                    for line in infile:
                        outfile.write(line)
        cmd_str = "{}/makeblastdb -in {} "\
                  "-parse_seqids -dbtype nucl".format(blast_bin, db_path)
        cmd = shlex.split(cmd_str)
        subprocess.call(cmd)
    print("Database is ready!")
    return db_path


def clean_up(temp_file):
    subprocess.call(["rm", temp_file])
    return os.path.isfile(temp_file)


def gene_nucleotide_blast(gene_fasta, blast_db, blast_bin, output_directory):
    """
    Running blast on a gene against a db using default params, except for evalue
    evalue: 1e-06
    :param gene_fasta:
    :param blast_db:
    :param blast_bin:
    :param output_directory:
    :return:
    """
    gene_name = os.path.basename(gene_fasta).split(".")[0]
    prefix = ut.to_str(gene_name) + "_blast.output.txt"
    outfile = os.path.join(output_directory, prefix)
    cmd_str = '{}/blastn -db {} -query {} '\
              '-evalue 0.000001 -outfmt '\
              '"6 qacc sacc qstart qend sstart send evalue pident qcovhsp"'.format(blast_bin,
                                                                                   blast_db,
                                                                                   gene_fasta)
    # IMPORTANT: if change output format need to change read_blast_output accordingly
    cmd = shlex.split(cmd_str)
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.read()
    with open(outfile, "w") as fo:
        fo.write(ut.to_str(output))
    return outfile


def get_sequence_out_of_multi_fasta(multi_fasta):
    """
    Generator yeilding one seq at a time from a multi_fasta file
    :param multi_fasta: path to multi_fasta file
    :yields: tuple of gene name, and gene sequence

    """
    with open(multi_fasta, "r") as fh:
        for line in multi_fasta:
            while len(line) > 0:
                gene_name = fh.readline()
                gene_seq = fh.readline()
                yield (gene_name, gene_seq)


def create_fasta_file(gene_name, gene_seq, output_directory, today):

    prefix = gene_name.strip().strip(">").replace(":", "_")
    temp_file = os.path.join(output_directory, ("{}_{}.fa".format(today, prefix)))
    with open(temp_file, "w") as tf:
        tf.write(gene_name)
        tf.write(gene_seq)
    return temp_file

def run_nucleotide_blast(output_directory, blast_bin, db_prefix, genome_folder,
                         multi_fasta, today=TODAY):

    # 1 Make/Check on database -> output db name
    db_path = make_nucleotide_blast_db(genome_folder, output_directory,
                                       db_prefix, blast_bin)

    # 2 Get gene from multifasta, save to file -> output file name
    multi_fasta_generator = get_sequence_out_of_multi_fasta(multi_fasta)
    gene_name, gene_seq = next(multi_fasta_generator)
    while gene_name:
        temp_file = create_fasta_file(gene_name, gene_seq, output_directory, today)
        # 3 Run blast on that file -> output blast results
        blast_out = gene_nucleotide_blast(temp_file, db_path, blast_bin, output_directory)
        # 4 Remove gene file
        if os.path.isfile(blast_out):
            clean_up(temp_file)
        gene_name, gene_seq = next(multi_fasta_generator)




def read_blast_output_file(blast_output_file):
    gene = ''
    blast_out = []
    gene_identity = []
    gene_coverage = []
    with open(blast_output_file, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                info = line.rstrip().split("\t")
                gene = info[0].split(":")[1]
                # this might change if using different assemblies
                contig = info[1]
                start = info[4]
                end = info[5]
                evalue = float(info[6])
                identity = float(info[7])
                coverage = float(info[8])
                if int(end) < int(start):
                    start, end = end, start
                if evalue < 0.000001 and identity > 80 and coverage > 90:
                    blast_out.append((contig, start, end))
                    gene_identity.append(identity)
                    gene_coverage.append(coverage)
        mean_identity = float(np.mean(gene_identity))
        mean_coverage = float(np.mean(gene_coverage))
    return gene, blast_out, round(mean_identity, 2), round(mean_coverage, 2)


def process_blast_output(blast_out_dir):

    """
    :param blast_out_dir:
    :return: location of blast hits within genomes

    """
    gene_in_genomes = {}
    gene_identities = {}
    gene_coverage = {}
    file_names = [os.path.join(blast_out_dir, f) for f in os.listdir(blast_out_dir)]
    for fi in file_names:
        if "blast.output" not in fi:
            continue
        gene, locations, mean_identity, mean_coverage = read_blast_output_file(fi)
        gene_in_genomes[gene] = locations
        gene_identities[gene] = mean_identity
        gene_coverage[gene] = mean_coverage
    return gene_in_genomes, gene_identities, gene_coverage


############################################################################################
# For each processed file need to find overlaps with gff

def find_overlap(gene_in_genome, gff):
    gff_bed = pybedtools.BedTool(gff)
    gene_loc_string = " ".join(gene_in_genome)
    gene_loc_bedtool = pybedtools.BedTool(gene_loc_string, from_string=True)
    intersect = gff_bed.intersect(gene_loc_bedtool, u=True)
    # u: Write original A entry once if any overlaps found in B.
    # In other words, just report the fact at least one overlap was found in B
    if not intersect:
        return "p{}_{}".format(gene_in_genome[1], gene_in_genome[2])
    df = intersect.to_dataframe()
    prokka = df.iloc[0].attributes
    prokka = prokka.split(";")[0].split("=")[1].strip()
    # Seems like a complicated way to get to prokka, is there a better way?
    return prokka


def find_all_overlaps(gene_in_genomes, gff_folder):
    gene_to_prokka = {}
    for gene, locations in gene_in_genomes.items():
        gene_to_prokka[gene] = {}
        for loc in locations:
            genome = loc[0].split("_")[0]  # This will have to change depending on assembly contig naming
            gff = os.path.join(gff_folder, "{}.gff".format(genome))  # This is specific to gff naming format
            if not os.path.isfile(gff):
                genome = "HM0"+genome.split("HM")[1]  # Awkward work around to avoid errors from HM1 vs HM01
                gff = os.path.join(gff_folder, "{}.gff".format(genome))
            prokka = find_overlap(loc, gff)
            gene_to_prokka[gene][genome] = prokka
    return gene_to_prokka



def make_presence_absence_matrix(gene_to_prokka, config_dict, out_dir,
                                 genome_set="all"):
    """

    :param gene_to_prokka: nested dics mapping gene to prokka, ex. {c333:{HM01:PROKKA1, HM03:PROKKA3}}
    :param config_dict:
    :param output_directory:
    :param genome_set: options right now just 'all'
    :return:
    """

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
    matrix_file = os.path.join(out_dir,
                               "{}_presence_absence_matrix.csv".format(TODAY))
    matrix_prokka_file = os.path.join(out_dir,
                                      "{}_presence_absence_matrix_prokka.csv".format(TODAY))
    matrix = pd.DataFrame(pa_matrix).T
    matrix_prokka = pd.DataFrame(pa_matrix_prokka).T
    matrix.to_csv(matrix_file)
    matrix_prokka.to_csv(matrix_prokka_file)
    return matrix_file


def get_gene_sequences_from_kegg(gene_list_file, genome="", out_dir=".", blast=""):
    """

    :param gene_list_file: txt file with locus tags, one on each line
    :param genome: KEGG genome designation, ex. ecc for CFT073
    :param output_directory: directory where output files will be saved
    :param blast: "nt" or "aa" whehter the function will return path to nucleotide or amino acid fasta,
    if nothing is specified will return path to the info file
    :return: one of three file names

    """
    if genome:
        gene_list = ["{}:{}".format(genome, gene.strip()) for gene in open(gene_list_file, "r")]
    else:
        gene_list = [gene.strip() for gene in open(gene_list_file, "r")]

    gene_set = GeneSet(gene_list, out_dir=out_dir)
    gene_set.get_info_for_each_gene()
    if blast == "nt":
        return gene_set.write_nt_seq()
    elif blast == "aa":
        return gene_set.write_aa_seq()
    else:
        return gene_set.get_info_df()



def compare_virulence_gene_expression(virulence_factor_file="", out_dir="",
                                      config="", blast="nt", multi_fasta="",
                                      blast_run=True, clean_up=False, today=TODAY):
    if not config:
        config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")
    config_dict = ut.process_config(config)
    if not out_dir:
        out_dir = config_dict["data"]["shared"]

    # 1 Get nucleotide sequence for each virulence factor
    if blast_run and not multi_fasta:
        if not virulence_factor_file:
            virulence_factor_file = config_dict["virulence_genes"]["path"]
        multi_fasta = get_gene_sequences_from_kegg(virulence_factor_file, '', out_dir, blast)


    # 2 For each gene run blast
    if blast_run:
        blast_bin = config_dict["blast"]["bin"]
        genome_dir = config_dict["blast"]["genome_dir"]
        db_prefix = config_dict["blast"]["db_prefix"]
        run_nucleotide_blast(out_dir, blast_bin, db_prefix, genome_dir, multi_fasta, today=TODAY)

    # 3 Process Blast output
    gene_in_genomes, gene_identities, gene_coverage = process_blast_output(out_dir)

    # 4 Identify Prokkas
    gff_folder = config_dict["blast"]["gff_dir"]
    gene_to_prokka = find_all_overlaps(gene_in_genomes, gff_folder)
    pd.DataFrame(gene_to_prokka).T.to_csv(os.path.join(out_dir,
                                          "{}_virulence_orthologs.csv".format(TODAY)))


    # 5 Build p/a matrix
    # make_presence_absence_matrix(gene_to_prokka, config_dict,# todo DONT think i need this
    #                              out_dir, genome_set="all")

    # 6 write out mean identity and coverage
    identity_file = os.path.join(out_dir, (today + "_identity_coverage.csv"))
    with open(identity_file, "w") as i_f:
        i_f.write('"",identity,coverage\n')
        for gene in gene_identities.keys():
            i_f.write('{},{},{}\n'.format(gene, gene_identities[gene], gene_coverage[gene]))
    # 7 # todo add cleanup step

if __name__ == "__main__":

    args = get_args().parse_args()
    if args.multifasta:
        mf = args.multifasta
    else:
        mf = ""
    if args.analysis == "blast" or args.analysis == "both":
        blast_run = True
    else:
        blast_run = False

    compare_virulence_gene_expression(virulence_factor_file=args.virulence_file,
                                      out_dir=args.outdir,
                                      config="", blast="nt", multi_fasta=mf,
                                      blast_run=blast_run, clean_up=False)


