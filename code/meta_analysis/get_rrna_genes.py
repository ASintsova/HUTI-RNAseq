import pandas as pd
import sys

sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/code/methods')
import keggAPI
import helpers


def get_rrna_genes(outdir, genome="eco", search_term="ribosomal RNA"):
    genes = keggAPI.searchKEGGGenome(genome, search_term, outdir)
    print(genes)


if __name__ == "__main__":
    odir = "/Users/annasintsova/git_repos/HUTI-RNAseq/results/meta_analysis"
    get_rrna_genes(odir)