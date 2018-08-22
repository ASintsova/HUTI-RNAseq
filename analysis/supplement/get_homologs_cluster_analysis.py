import os
import subprocess
import shlex
from collections import Counter
import re
import glob
import pandas as pd
import itertools
import shutil
import sys
sys.path.append("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/methods")
import run_gethomologues as gh

def cluster_size_distribution(cluster_list):
    cluster_sizes = []
    with open(cluster_list, "r") as cl:
        for line in cl:
            result = re.search("size=\d*", line)
            if result:
                cluster_sizes.append(int(result.group().split("=")[1]))
    count_clusters = Counter(cluster_sizes)
    return count_clusters


def cluster_size_distribution_intersection(pangenome_matrix):
    # ???
    df = pd.read_csv(pangenome_matrix, sep="\t", index_col=0)
    return Counter(list(df.sum(0)))


def cluster_sizes(run_id):
    # run_id: path to the output directory
    cog_cluster_list = glob.glob("genomes_homologues/*algCOG*.cluster_list")[0]
    omcl_cluster_list = glob.glob("genomes_homologues/*algOMCL*.cluster_list")[0]
    pangenome_matrix = glob.glob("genomes_homologues/*pangenome_matrix*.tab")[0]
    cog_clusters = cluster_size_distribution(cog_cluster_list)
    omcl_clusters = cluster_size_distribution(omcl_cluster_list)
    intersection_clusters = cluster_size_distribution_intersection(pangenome_matrix)
    clusters = [cog_clusters, omcl_clusters, intersection_clusters]
    df = pd.concat([pd.DataFrame.from_dict(d, orient='index') for d in clusters], axis=1)
    df.columns = ["COG", "OMCL", "Intersection"]
    df.to_csv("{}_cluster_sizes.csv".format(run_id))


def params(C = [50, 60, 70, 75] , S =[1, 70, 90] ):
    return list(itertools.product(C, S))


def generate_cluster_size_distributions(gh, gbk_dir, I, params, exclude_paralogs):

    C = params[0]
    S = params[1]
    e_name = "e1" if exclude_paralogs else "e0"

    run_id = "C{}_S{}_{}_".format(C, S, e_name)
    gh.run_get_homologs(gh, gbk_dir, I, C, S, exclude_paralogs, run_id)
    cluster_sizes(run_id)


if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("lib/config")
    gh = config.get("GETHOMS", "bin")
    gbk_dir = config.get("GETHOMS", "gbk_dir")
    I = config.get("GETHOMS", "I")

    # params = params()
    # C = [50, 60, 70, 75], S = [1, 70, 90]
    params = (75, 1)
    # generateClusterSizeDistributions(gh, gbk_dir, I, params, exclude_paralogs)
    ###### Include Paralogs ######

    ###### Only Syntenic genes (inlcuding paralogs)######