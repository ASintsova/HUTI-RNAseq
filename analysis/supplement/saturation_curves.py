import numpy as np
import os
import pandas as pd
import pybedtools
import pysam
import random
import sys

sys.path.append("/Users/annasintsova/tools/bedtools2/bin/bedtools")
sys.path.append("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/methods")

import helpers


def subsample_bam(bam_file, subsampled_bam, fraction):
    """
    Subsample without replacement using pysam

    """

    bam = pysam.AlignmentFile(bam_file)
    output = pysam.AlignmentFile(subsampled_bam, "wb", template=bam)  # subsampled file
    for read in bam.fetch():
        if random.random() < fraction:
            output.write(read)
    bam.close()
    print("{} subsampled".format(bam_file))
    output.close()


def coverage_pybedtools(bam_file, gff_file, coverage_file):
    gff = pybedtools.BedTool(gff_file)
    bam = pybedtools.BedTool(bam_file)
    coverage = gff.coverage(bam)
    coverage.saveas(coverage_file)


def high_covered_genes(coverage_file, cutoff, reads=True):

    """
    Resulting coverage file columns: first 9 are the gff columns,
        10: number of reads overalapped with gene,
        11: number of bases in the gene covered by reads
        12: length of gene
        13: % of gene covered

    """

    df = pd.read_csv(coverage_file,  sep="\t", header=None,
                     names=["chr", "soft", "feat", "start", "end",
                            "score", "str", "frame", "attr", "num_reads",
                            "num_bases_covered", "len_gene", "coverage"])
    if reads:
        return df.loc[df["num_reads"] >= cutoff].shape[0]
    else:
        return df.loc[df["coverage"] >= cutoff].shape[0]


def find_gene_num_in_iter(bam_file, subsample_path, gff_file,
                          coverage_path, fraction, cutoff, reads,  n):
    """
    Runs through subsampling and coverage calculation n times, returns list of num of hi genes
    for each iteration

    """
    def num_genes_in_one_iter():
        for i in range(n):
            subsample_bam(bam_file, subsample_path, fraction)
            coverage_pybedtools(subsample_path, gff_file, coverage_path)
            num_high_genes = high_covered_genes(coverage_path, cutoff, reads)
            yield num_high_genes
    return np.mean(list(num_genes_in_one_iter()))


def calculate_saturation_curve_data(bam_file, gff_file, subsample_range,
                                    cutoff, reads, output_dir, n):

    def generate_saturation_curve_data_for_fraction():
        for fraction in subsample_range:
            prefix = helpers.to_str(os.path.basename(bam_file).split(".")[0])
            subsampled_bam = prefix + "_{}_subsample.bam".format(str(fraction))
            subsample_path = os.path.join(output_dir, subsampled_bam)
            coverage_file = prefix + "_{}_coverage.txt".format(str(fraction))
            coverage_path = os.path.join(output_dir, coverage_file)
            mean_num_genes = find_gene_num_in_iter(bam_file, subsample_path, gff_file,
                                                   coverage_path, fraction, cutoff, reads, n)
            os.remove(subsample_path)
            os.remove(coverage_path)
            yield (fraction, mean_num_genes)
    return list(generate_saturation_curve_data_for_fraction())


def get_samples(config_dict, strains='all', conditions='all'):
    strains = config_dict["strains"][strains].split()
    conditions = config_dict["conditions"][conditions].split()
    samples = ["{}_{}".format(s, t) for s in strains for t in conditions]
    return samples


def saturation_curves(config, subsample_range, cutoff, reads, i, strains='all', conditions='all'):
    config_dict = helpers.process_config(config)
    samples = get_samples(config_dict, strains, conditions)
    out_dir = config_dict["output_directory"]["results"]
    df = pd.DataFrame(index=[str(c * 100) + "%" for c in subsample_range])
    for sample in samples:
        genome = sample.split("_", 1)[0]
        bam_path = config_dict["BAMS"]["path"]
        bam_file = os.path.join(bam_path, "{}_trimmed_sorted.bam".format(sample))
        gff_path = config_dict["GFF"]["path"]
        gff_file = os.path.join(gff_path, "{}.gff".format(genome))
        output_file = os.path.join(out_dir, "{}_saturation_curve_data.csv".format(sample))
        results = calculate_saturation_curve_data(bam_file, gff_file, subsample_range,
                                                  cutoff, reads, out_dir, i)
        df[sample] = [y[1] for y in results]
        df.to_csv(os.path.join(output_file))


if __name__ == "__main__":
        config_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/supplement/config"
        fractions = [0.001, 0.005, 0.01, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        hi_reads = 10
        r = True
        iterations = 5
        saturation_curves(config_file, fractions, hi_reads, r, iterations, 'all', 'all')

        # Same problem as with gh program - need it to spit out intermediate results
        # Also very slow :'(
