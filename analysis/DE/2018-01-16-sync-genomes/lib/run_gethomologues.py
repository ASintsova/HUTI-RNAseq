"""


parameters I want to be able to change:

-C 50/60/70/75(default)
-S  1(default)/70/90
-e exclude paralogs

always to include
 -A
 -c


In the end want to check for # of singletons: go through all files in intersect_pan_... and count how many files have only one entry, do a counting dictionary can then see distribution



Do I want to do this on smaller number of genomes? Yes, let's do 5: CFT073, HM01, HM07, HM57, HM66

Do I want to adapt this to cluster scenario? Not yet.


"""

import configparser
import os
import subprocess
import shlex
from collections import Counter
import re
import glob
import pandas as pd
import itertools
import shutil
import threading

#The gh output directory will be the same as the one from which this script is run

"""
        for compare clusters also want to try -s (synteny) option and maybe -I (produce clusters with single-copy seqs from ALL taxa in file) -m is for pangenomic matrices and -T for the tree
        which I actually don't really need that
        
"""

def clusterSizeDistribution(cluster_list):
    cluster_sizes = []
    with open(cluster_list, "r") as cl:
        for line in cl:
            result = re.search("size=\d*", line)
            if result:
                cluster_sizes.append(int(result.group().split("=")[1]))
    count_clusters = Counter(cluster_sizes)
    return count_clusters

def clusterSizeDistributionIntersection(pangenome_matrix):
    df = pd.read_csv(pangenome_matrix, sep="\t", index_col=0)
    return Counter(list(df.sum(0)))


def clusterSizes(run_id): #path to the output directory

    algCOG_cluster_list = glob.glob("genomes_homologues/*algCOG*.cluster_list")[0]
    algOMCL_cluster_list = glob.glob("genomes_homologues/*algOMCL*.cluster_list")[0]
    pangenome_matrix = glob.glob("genomes_homologues/*pangenome_matrix*.tab")[0]

    cog_clusters = clusterSizeDistribution(algCOG_cluster_list)
    omcl_clusters = clusterSizeDistribution(algOMCL_cluster_list)
    intersection_clusters = clusterSizeDistributionIntersection(pangenome_matrix)

    df = pd.concat([pd.DataFrame.from_dict(d, orient='index') for d in [cog_clusters, omcl_clusters, intersection_clusters]], axis=1)
    df.columns = ["COG", "OMCL", "Intersection"]
    df.to_csv("{}cluster_sizes.csv".format(run_id))


def runGH(gh, gbk_dir, I, C, S, exclude_paralogs, run_id):
    # Make a temp directory and go there
    if not os.path.exists(run_id):
        os.makedirs(run_id)
    os.chdir(run_id)
    # All other paths should be absolute, including I

    # Copy gbk files specified in I into a new genome directory
    os.makedirs('genomes')
    with open(I, "r") as gf:
        for line in gf:
            print(line)
            file = os.path.join(gbk_dir, line.strip())
            shutil.copy(file, 'genomes')

    #Run gh with appropriate params
    ####################################################

    #Now there's not I, and gbk dir is hard coded
    e = " -e" if exclude_paralogs else ''
    e_name = "1" if exclude_paralogs else "0"

    cmd1_str = "{} -d genomes -G -t 0 -C {} -S {} -A -c{}".format(os.path.join(gh, "get_homologues.pl"), C, S, e) # Running COGtriangle algorithm
    print(cmd1_str)
    cmd1 = shlex.split(cmd1_str)
    subprocess.call(cmd1)
    cmd2_str = "{} -d genomes -M -t 0 -C {} -S {} -A -c{}".format(os.path.join(gh, "get_homologues.pl"), C, S, e)
    cmd2 = shlex.split(cmd2_str)
    subprocess.call(cmd2)
    gh_out_dir = "genomes_homologues"

    cog_out = glob.glob(gh_out_dir+"/*algCOG_e{}_C{}_S{}_".format(e_name, str(C), str(S)))[0]
    omcl_out = glob.glob(gh_out_dir+"/*algOMCL_e{}_C{}_S{}_".format(e_name, str(C), str(S)))[0]

    if cog_out.split("alg")[0] == omcl_out.split("alg")[0]:
        print("Go!")
        #s = " -s" if synteny else ''
        cmd3_str = "{} -o {}/run_{}_pan_C{}_S{} -d {},{},  -t 0 -m -T".format(os.path.join(gh, "compare_clusters.pl"), gh_out_dir, run_id, str(C), str(S), cog_out, omcl_out)
        print(cmd3_str)
        cmd3 = shlex.split(cmd3_str)
        print(cmd3)
        subprocess.call(cmd3)

    ####################################################
    #Delete everything except for .tab and _list files
    for root, dirs, files in os.walk("genomes_homologues", topdown=False):
        for name in files:
            if '.tab' in name or '_list' in name:
                print(os.path.join(root, name))
            else:
                os.remove(os.path.join(root, name))
        for name in dirs:
            fold = os.path.join(root, name)
            print(fold)
            for f in os.listdir(fold):
                shutil.move(os.path.join(fold, f), root)

        for name in dirs:
            if not os.listdir(os.path.join(root, name)):
                os.removedirs(os.path.join(root, name))


def params(C = [50, 60, 70, 75] , S =[1, 70, 90] ):
    return list(itertools.product(C, S))


def generateClusterSizeDistributions(gh, gbk_dir, I, params, exclude_paralogs):

    C = params[0]
    S = params[1]
    e_name = "e1" if exclude_paralogs else "e0"

    run_id = "C{}_S{}_{}_".format(C, S, e_name)
    runGH(gh, gbk_dir, I, C, S, exclude_paralogs, run_id)
    clusterSizes(run_id)

# def generateClusterSizeDistributions(params, I, exclude_paralogs, gh, gbk_dir):
#     e_name = "1" if exclude_paralogs else "0"
#     gh_dir = os.path.basename(gbk_dir)+ "_homologues"
#     results = {}
#     print(os.path.isdir(gh_dir))
#     for i in range(len(params)):
#         #run_id = "id" + str(i)
#
#         C = params[i][0]
#         S = params[i][1]
#         run_id = "id_C{}_S{}".format(str(C), str(S))
#         runGH(gh, gbk_dir, I, C, S, exclude_paralogs, run_id)
#         COG_cl = glob.glob(gh_dir+"/*algCOG_e{}_C{}_S{}_.cluster_list".format(e_name, str(C), str(S)))[0]
#
#         cog_sizes = clusterSizeDistribution(COG_cl)
#         OMCL_cl = glob.glob(gh_dir+"/*algOMCL_e{}_C{}_S{}_.cluster_list".format(e_name, str(C), str(S)))[0]
#
#         #need to add checks and balances here
#         ocml_sizes = clusterSizeDistribution(OMCL_cl)
#         pangenome_matrix = glob.glob(gh_dir+"/run_{}_pan_C{}_S{}/pangenome_matrix*.tab".format(run_id, str(C), str(S)))[0]
#
#         pangenome_sizes = clusterSizeDistributionIntersection(pangenome_matrix)
#         #results[run_id] = {"COG": cog_sizes, "OCML": ocml_sizes, "pan": pangenome_sizes,
#         #                   "conditions": "C:{}, S: {}".format(str(C), str(S))}
#         results[run_id] = {"COG": cog_sizes, "OCML": ocml_sizes, "pan": pangenome_sizes}
#         #with open ("data/{}_intermediate.csv")
#         subprocess.call(["rm", "-r", gh_dir])
#
#     dd = {k: pd.DataFrame(v) for k, v in results.items()}
#     df = pd.concat(dd, axis=0)
#     df.to_csv("data/get_homologues_cluster_size_analysis.csv")


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read("lib/config")
    gh = config.get("GETHOMS", "bin")
    gbk_dir = config.get("GETHOMS", "gbk_dir")
    I = config.get("GETHOMS", "I")

    #params = params()
    #C = [50, 60, 70, 75], S = [1, 70, 90]
    params = (75, 1)
    #generateClusterSizeDistributions(gh, gbk_dir, I, params, exclude_paralogs)
    ###### Include Paralogs ######

    ###### Only Syntenic genes (inlcuding paralogs)######





"""

a few problems with this that need to be fixed:
way too slow
need to print out where we are in the run
need to save intermediate results, because I have no patience to wait for it to be done
re-running blast each time is very inefficient
"""