
import os
import pandas as pd

def full_crossRef(crossRef="/Users/annasintsova/git_repos/HUTI-RNAseq/"\
                        "data/get_homologs_output/C50_S90_e0_"\
                        "/run_C50_S90_e0__pan_C50_S90/2018-02"\
                            "-26_pangenome_matrix_t0_crossRef.csv"):

    with open(crossRef, "r") as fh:
        keys = fh.readline().split(",")[1:]
        cross_dict = {k:[] for k in keys}
        for line in fh:
            vals = line.split(',')[1:]
            for k, v in zip(keys, vals):
                cross_dict[k] += [v]

    return cross_dict


def get_homologues(genome,
                   cross_dict):
    """
    :param genome: reference genome
    :param crossRef: crossRef file
    :return: {Prokka1:[(HM3, Prokka56), (HM54, Prokka65)]}
    """

    genes = cross_dict[genome]
    keys = [k for k in cross_dict.keys() if "_" not in k and k != genome]
    print(keys)
    genome_homologues = {}
    for i in range(3302, 3308):
        if genes[i]:
            print(genes[i])
            genome_homologues[genes[i]] = [(k, cross_dict[k][i]) for k in keys if cross_dict[k][i]]
    print(genome_homologues)


def get_counts(genome, cross_dict):
    print("Get counts")


def get_RPKMS(genome, rpkms):
    print("Get RPKMS")

if __name__ == "__main__":
    cross_dict = full_crossRef()
    print(get_homologues("HM54", cross_dict))