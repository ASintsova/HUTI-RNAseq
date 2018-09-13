# todo not sure this actually works


import sys
sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/')
sys.path.insert(0, "/Users/annasintsova/tools/bedtools2/bin/")
from methods import saturation_curves



####################
# Running through saturation_curves functions

def test_saturation_curves(bam_file, subsample_bam, gff_file, fraction, coverage_file, coverage_file2):
    saturation_curves.subSampleBam(bam_file, subsample_bam, fraction)
    saturation_curves.coverageBedTools(subsample_bam, gff_file, coverage_file)
    saturation_curves.coverageBedTools(subsample_bam, gff_file, coverage_file)
    saturation_curves.coverage_pybedtools(subsample_bam, gff_file, coverage_file2)

if __name__ == "__main__":
    bam_file = "code/tests/alignments/HM01_UTI_trimmed_sorted.bam"
    subsample_bam = "code/tests/output/subsample_test.bam"
    gff_file = "code/tests/annotations/HM1_gh_final.gff"
    coverage_file = "code/tests/output/coverage_test.bed"
    coverage_file2 = "code/tests/output/coverage_test2.bed"
    fraction = 0.05
    test_saturation_curves(bam_file, subsample_bam,
                           gff_file, fraction, coverage_file,
                           coverage_file2)





