import numpy as np
import sys
sys.path.append("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/supplement")

import saturation_curves as sc

def test_find_gene_num_in_iter():
    #bam_file, gff_file = bam_gff
    bam_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/alignments/HM66_UTI_trimmed_sorted.bam"
    gff_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/gff_files/HM66.gff"
    subsample_path = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/test_data/test.bam"
    coverage_path = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/test_data/test.bed"
    fraction = 0.05
    cutoff = 10
    reads = True
    n = 5
    nums = sc.find_gene_num_in_iter(bam_file, subsample_path, gff_file,
                                    coverage_path, fraction, cutoff, reads, n)
    print(np.mean(nums), np.std(nums))


def test_calculate_saturation_curve_data():
    bam_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/alignments/HM66_UTI_trimmed_sorted.bam"
    gff_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/gff_files/HM66.gff"
    output_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/test_data/"
    subsample_range = [0.05, 0.1]
    cutoff = 10
    reads = True
    n = 5
    results = sc.calculate_saturation_curve_data(bam_file, gff_file, subsample_range,
                                       cutoff, reads, output_dir, n)
    print(results)

def test_get_samples():
    config = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/supplement/config"
    s = sc.get_samples(config, strains='all', conditions='all')
    print(s)
if __name__ == "__main__":
    test_get_samples()