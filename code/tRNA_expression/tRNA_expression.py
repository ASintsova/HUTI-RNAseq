from collections import namedtuple
import datetime as dt
import os
import shlex
import subprocess


def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode("utf-8")
    else:
        value = bytes_or_str
    return value


def find_trnas(fasta_folder, output_directory):

    """
    find tRNAs!
    :return:
    """
    files = [os.path.join(fasta_folder, c) for c in os.listdir(fasta_folder)]
    for fasta_file in files:
        suffix = os.path.basename(fasta_file).split('.')[0] +"_.tRNAscan" # todo import helpers, convert to str
        output_file = os.path.join(output_directory, suffix)
        cmd_str = "/Users/annasintsova/Downloads/tRNAscan-SE-2.0/tRNAscan-SE" \
                  " -B {}".format(fasta_file)
        cmd = shlex.split(cmd_str)
        print(cmd)
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        output = output.stdout.read()

        with open(output_file, "w") as fo:
            fo.write(to_str(output))

def proccess_trna_scan_output(trna_output_file, output_file):

    trna = namedtuple('trna', 'seq start end id')
    all_trnas = []
    trnas_count = {}

    with open(trna_output_file, "r") as fh:
        for line in fh:
            if not line.startswith("HM"):
                continue
            else:
                words = line.strip().split('\t')
                seq = words[0].strip()
                loc1 = int(words[2])
                loc2 = int(words[3])
                species = words[4] + "_" + words[5]
                trnas_count[species] = trnas_count.get(species, 0) + 1
                print(trnas_count)
                name = species +"_" + str(trnas_count[species])
                all_trnas.append(trna(seq, min(loc1, loc2), max(loc1, loc2), name))
    with open(output_file, "w") as fo:
        for tr in all_trnas:
            print(tr)
            fo.write("{}\t{}\t{}\tid={}\n".format(tr.seq, tr.start, tr.end, tr.id))


def find_expression_with_bedtools(gff, bam ):
    cmd_str = "bedtools coverage -counts -a {} -b {}".format(gff, bam)

    cmd = shlex.split(cmd_str)
    print(cmd)
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = output.stdout.read()
    output_file = gff.split('.gff')[0] + ".counts"
    with open(output_file, "w") as fo:
        fo.write(to_str(output))

if __name__ == "__main__":
    trna_output_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/code/tRNA_expression/tests/test_data/HM01_.tRNAscan"
    gff = '/Users/annasintsova/git_repos/HUTI-RNAseq/code/tRNA_expression/tests/test_data/HM01.gff'
    bam = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/alignments/HM01_UTI_trimmed_sorted.bam"
    bam2 ="/Users/annasintsova/git_repos/HUTI-RNAseq/data/alignments/HM01_UR_trimmed_sorted.bam"
    #proccess_trna_scan_output(trna_output_file, gff)
    find_expression_with_bedtools(gff, bam)