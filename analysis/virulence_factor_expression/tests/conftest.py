import datetime as dt
import os
import pytest
import sys
import shutil

sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/methods')
import helpers


@pytest.fixture()
def setup(tmpdir):
    config_dict = helpers.process_config("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis"
                                "/virulence_factor_expression/config")

    virulence_file = config_dict["virulence_genes"]["path"]
    output_directory = str(tmpdir)
    yield virulence_file, output_directory
