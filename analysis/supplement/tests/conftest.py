import pytest

@pytest.fixture()
def bam_gff():
    bam_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/alignments/HM66_UTI_trimmed_sorted.bam"
    gff_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/gff_files/HM66.gff"

    yield bam_file, gff_file
