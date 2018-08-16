import pybedtools
import datetime as dt
import os
import shlex
import subprocess
import sys

sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/code/methods')

import helpers

TODAY = dt.datetime.today().strftime("%Y-%m-%d")

def process_config():
    print("Find similar function elsewhere")


##################################################################################
# RUN NUCLEOTIDE BLAST

def make_nucleotide_blast_db(genome_folder, output_directory, db_prefix,
                             blast_bin):

    """

    :param genome_folder:
    :param output_directory:
    :param db_prefix:
    :param blast_bin:
    :return:

    """
    # Check if database exists in output directory
    genomes = [os.path.join(genome_folder, g) for g in os.listdir(genome_folder)]
    db_path = os.path.join(output_directory, (db_prefix +".fsa"))
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
    prefix = os.path.basename(gene_fasta).split(".")[0] + "_blast.output.txt"
    outfile = os.path.join(output_directory, prefix)
    cmd_str = '{}/blastn -db {} -query {} -outfmt '\
               '"7 qacc sacc evalue qstart qend sstart send nident"'.format(blast_bin,
                                                                     blast_db,
                                                                     gene_fasta)
    cmd = shlex.split(cmd_str)
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.read()
    with open(outfile, "w") as fo:
        fo.write(helpers.to_str(output))
    return outfile


def get_sequence_out_of_multi_fasta(multi_fasta=''):

    with open(multi_fasta, "r") as fh:
        for line in multi_fasta:
            while len(line) > 0:
                gene_name = fh.readline()
                gene_seq = fh.readline()
                yield (gene_name, gene_seq)


def run_nucleotide_blast(output_directory, blast_bin, db_prefix, genome_folder,
                         multi_fasta, today=TODAY):

    #1 Make/Check on database -> output db name
    db_path = make_nucleotide_blast_db(genome_folder, output_directory, db_prefix,
                             blast_bin)
    #2 Get gene from multifasta, save to file -> output file name
    multi_fasta_generator = get_sequence_out_of_multi_fasta(multi_fasta)
    gene_name, gene_seq = next(multi_fasta_generator)
    while gene_name:
        # 2 Get gene from multifasta, save to file -> output file name
        prefix = gene_name.strip().strip(">").replace(":", "_")
        temp_file = os.path.join(output_directory,
                                 ("{}_{}.fa".format(today, prefix)))
        with open(temp_file, "w") as tf:
            tf.write(gene_name)
            tf.write(gene_seq)
        #3 Run blast on that file -> output blast results
        blast_out = gene_nucleotide_blast(temp_file, db_path, blast_bin, output_directory)
        if os.path.isfile(blast_out):
        #4 Remove gene file
            clean_up(temp_file)
        gene_name, gene_seq = next(multi_fasta_generator)

def clean_up(temp_file):

    subprocess.call(["rm", temp_file])
    return os.path.isfile(temp_file)


###################################################################################


#For each gene there was a blast.out generated, read it in and delete it

def process_blast_output(blast_out_dir):

    """
    ONLY LOOKING AT FULL LENGTH ALIGNMENTS
    :param blast_out_dir:
    :return: location of blast hits within genomes


    """
    gene_in_genomes = {}
    genomes_with_genes = {}
    file_names = [os.path.join(blast_out_dir, f) for f in os.listdir(blast_out_dir)]

    for fi in file_names:
        if not "blast.output" in fi:
            continue
        with open(fi, "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                else:
                    info = line.rstrip().split("\t")
                    gene = info[0].split(":")[1]
                    genome = info[1].split("_")[0]
                    sequence = info[1]
                    end = info[-1]
                    start = info[-2]

                    if info[3].strip() != '1':
                        continue
                    if int(end) < int(start):
                        start, end = end, start

                    if not gene in gene_in_genomes.keys():
                        gene_in_genomes[gene] = [(sequence,start, end)]
                    else:
                        gene_in_genomes[gene] += [(sequence, start, end)]
    return gene_in_genomes


####################################################################################
# for each processed file need to find overlaps with gff

def find_overlap(gene_in_genome, gff):
    a = pybedtools.BedTool(gff)
    b_string = " ".join(gene_in_genome)
    b = pybedtools.BedTool(b_string, from_string=True)
    intersect = a.intersect(b, u=True)
    if not intersect:
        return "p{}_{}".format(gene_in_genome[1], gene_in_genome[2])
    df = intersect.to_dataframe()
    prokka = df.iloc[0].attributes
    prokka = prokka.split(";")[0].split("=")[1].strip()
    return prokka

def find_all_overlaps(gene_in_genomes, gff_folder):
    gene_to_prokka = {}

    for gene, info in gene_in_genomes.items():
        gene_to_prokka[gene] = []
        for gene_in_genome in info:
            genome = gene_in_genome[0].split("_")[0]
            gff = os.path.join(gff_folder, "{}.gff".format(genome))
            if not os.path.isfile(gff):
                genome = "HM0"+genome.split("HM")[1]
                gff = os.path.join(gff_folder, "{}.gff".format(genome))
            prokka = find_overlap(gene_in_genome, gff)
            gene_to_prokka[gene] += [(genome, prokka)]

    return gene_to_prokka


# def extract sequence
# def align sequences



if __name__ == "__main__":

    output_directory = "/Users/annasintsova/git_repos/HUTI-RNAseq/code/tests/output"
    genome = "ecc"
    gene_list_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/code/tests/gene_list.txt"
    genome_folder = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/genomes/clinical_strains"
    blast_bin = "/Users/annasintsova/tools/ncbi-blast-2.7.1+/bin/"
    gene_fasta = "/Users/annasintsova/git_repos/HUTI-RNAseq/code/tests/output/gene_list_nt.fasta"
    #get_gene_blast_results(gene_list_file, genome,output_directory)
    #db_path = make_nucleotide_blast_db(genome_folder, output_directory, "test_db", blast_bin)
    #run_nucleotide_blast(gene_fasta, db_path, blast_bin, output_directory)
    #run_nucleotide_blast(output_directory, blast_bin, "test_db", genome_folder,
     #                   gene_fasta, today=TODAY)
    x = process_blast_output(output_directory)

    gff = "/Users/annasintsova/Downloads/RNA_reference_genomes/"
    print(find_all_overlaps(x, gff))



