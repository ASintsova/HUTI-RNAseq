import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
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


def get_labels(df):
    """ df.index needs to be in form strain_condition"""
    labels = []
    for i in df.index:
        label = sample_label_from_sample_name(i)
        labels.append(label)
    return labels


def sample_label(strain, treat, jn=" | "):
    return strain + jn + treat + " "


def sample_label_from_sample_name(sample_name):
    strain = sample_name.split("_")[0]
    condition = ur if sample_name.split("_")[1] == "UR" else uti
    return sample_label(strain, condition)


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



"""
Plotting PCA elipses:

__author__:

"""

def plot_point_cov(points, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, nstd, ax, **kwargs)


def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip
