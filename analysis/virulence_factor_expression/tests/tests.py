import sys
import os
sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/methods')
sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/virulence_factor_expression')
import helpers
import virulence_factor_expression


def test_generate_multifasta(setup):
    file, out_dir = setup
    multi_fasta_file = virulence_factor_expression.generate_multifasta(file, out_dir, "nt")
    assert os.path.isfile(multi_fasta_file)


