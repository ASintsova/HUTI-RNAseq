import datetime as dt
import os
import pytest
import sys
import shutil
sys.path.append('.')
import workflow
from modules import helpers
@pytest.fixture()
def locus_tag():
    return "b0356"

@pytest.fixture()
def gene_name():
    return "frmA"

