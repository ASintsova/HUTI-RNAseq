import pybedtools
import datetime as dt
import numpy as np
import os
import shlex
import subprocess
import sys

sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/methods')

import helpers

TODAY = dt.datetime.today().strftime("%Y-%m-%d")


##################################################################################
# RUN NUCLEOTIDE BLAST

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
    prefix = helpers.to_str(gene_name) + "_blast.output.txt"
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
        fo.write(helpers.to_str(output))
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


def clean_up(temp_file):
    subprocess.call(["rm", temp_file])
    return os.path.isfile(temp_file)


###################################################################################


# For each gene there was a blast.out generated, read it in and delete it


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
                eval = float(info[6])
                identity = float(info[7])
                coverage = float(info[8])
                if int(end) < int(start):
                    start, end = end, start
                if eval < 0.000001 and identity > 80 and coverage > 90:
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
        if not "blast.output" in fi:
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
        gene_to_prokka[gene] = []
        for loc in locations:
            genome = loc[0].split("_")[0]  # This will have to change depending on assembly contig naming
            gff = os.path.join(gff_folder, "{}.gff".format(genome))  # This is specific to gff naming format
            if not os.path.isfile(gff):
                genome = "HM0"+genome.split("HM")[1]  # Awkward work around to avoid errors from HM1 vs HM01
                gff = os.path.join(gff_folder, "{}.gff".format(genome))
            prokka = find_overlap(loc, gff)
            gene_to_prokka[gene] += [(genome, prokka)]
    return gene_to_prokka


# def extract sequence
# def align sequences


if __name__ == "__main__":

    od = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/test_data/temp"
    genome_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/" \
                 "genomes/clinical_strains"
    b_db = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/virulence_factor_expression"
    b_bin = "/Users/annasintsova/tools/ncbi-blast-2.7.1+/bin/"
    mf = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/virulence_factor_expression/" \
         "VIRULENCE_FACTORS_nt.fasta"
    gff = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/gff_files"
    run_nucleotide_blast(od, b_bin, b_db, genome_dir, mf)
    q, x, y = process_blast_output(od)
    gp = find_all_overlaps(q, gff)
    print(gp)






