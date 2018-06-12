import os
import pybedtools
import sys
import time
sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/code/methods')
import helpers

def calculate_gene_codon_usage(fasta, codon, abundance):
    """

    :param fasta: gene sequence
    :param codon: codon of interest
    :param abundance: mRNA abundance
    :return: Cij * Ejk
    """
    cnt = 0
    assert len(fasta) % 3 == 0
    print(len(fasta))
    for i in range(0, len(fasta), 3):
        if fasta[i:i+3] == codon:
            cnt += 1
    return cnt*abundance


def yieled_fasta_from_gff(fasta, gff):
    """
    Given a genome and a gff file, yields one fasta sequence at a time
    Does not look at strandiness, etc
    :param fasta:
    :param gff:
    :return:
    """
    a = pybedtools.BedTool(gff)
    a = a.sequence(fi=fasta)
    with open(a.seqfn) as fn:
        for line in fn:
            if line.startswith(">"):
                continue
            else:
                yield line

def reverse_trancscribe(gene_seq):
    d = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rev_tran = [d[N] for N in gene_seq[::-1]]
    return "".join(rev_tran)


def get_transcript_from_gff(gff, fasta):

    a = pybedtools.BedTool(gff)
    seq = yieled_fasta_from_gff(fasta, gff)
    for feature in a.features():
        gene_seq = next(seq).strip()
        if feature.strand == '-':
            gene_seq = reverse_trancscribe(gene_seq)

        yield (feature.name, gene_seq.replace("T", "U"))

def codon_frequency_in_genome(codon, gff, fasta):
    total_codons = 0
    codon_cnts = 0
    for name, seq in get_transcript_from_gff(gff, fasta):
        print(name)
        total_codons += len(seq)/3
        for i in range(0, len(seq), 3):
            if seq[i:i+3] == codon:
                codon_cnts += 1
    print(total_codons)
    print(codon_cnts)
    return codon_cnts/total_codons*100

def get_counts(sample, config_dict):
    count_dict = {}

    cnt_path = config_dict["counts"]["path"]
    cnts_file = os.path.join(cnt_path, sample+".csv")
    with open(cnts_file, "r") as fh:
        fh.readline()
        for line in fh:

            words = line.strip().split(",")
            count_dict[words[0]] = float(words[1])
    return count_dict



def codon_frequency_in_transcriptome(codon, sample, config_dict):

    genome = sample.split("_")[0]
    gff_suffix = genome + ".gff"
    fasta_suffix = genome + ".fasta"

    count_dict = get_counts(sample, config_dict)
    gff = os.path.join(config_dict["gff_folder"]["path"], gff_suffix)
    fasta = os.path.join(config_dict["genomes"]["path"], fasta_suffix)

    total_codons = 0
    codon_cnts = 0
    missing = 0
    E = 0
    for name, seq in get_transcript_from_gff(gff, fasta):
        codon_in_gene = 0

        total_codons += len(seq)/3
        for i in range(0, len(seq), 3):
            if seq[i:i+3] == codon:
                codon_cnts += 1
                codon_in_gene += 1
        if name in count_dict.keys():
            cnts = count_dict[name] + 4
        else:
            cnts = 4
            missing +=1
        #print("codon in gene ", codon_in_gene)

        E += cnts*codon_in_gene
        #print("E",E)


    print(total_codons)
    print(codon_cnts)
    print(codon_cnts/total_codons)
    print("Missing ", missing)
    return E



if __name__ == "__main__":
    config_dict = helpers.process_config("config")
    print(codon_frequency_in_transcriptome("CCU", "HM56_UTI", config_dict))
