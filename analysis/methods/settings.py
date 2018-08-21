import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from operator import itemgetter
import os
import pandas as pd
import scipy
from scipy import stats
import seaborn as sns;sns.set_style("white")

today = dt.datetime.today().strftime("%Y-%m-%d")


def invnorm(x):
    return stats.norm.ppf((x.rank() - 0.5)/x.count())


strain_info_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/huti_patient_info_short.csv"


def sample_label(strain, treat, jn=" | "):
    return strain + jn + treat


col1 = "#f34236"
col2 = "#d6c571"
col3 = "#88bc67"
col4 = "#2e8174"
col5 = "#143969"

clrs = [col1, col2, col3, col4, col5]

my_cmap = LinearSegmentedColormap.from_list('custom blue', [col5, col4, col2, col1], N=256)

ur = "URINE"
uti = "PATIENT"
join = " | "

strain_qual = {'good': ["HM56", "HM14", "HM43", "HM54", "HM86"],
               'okay': ["HM56", "HM14", "HM43", "HM54", "HM86", "HM01", "HM03", "HM06", "HM68"],
               'so-so': ["HM57", "HM17", "HM07" "HM60"],
               'bad': ["HM66"]}
