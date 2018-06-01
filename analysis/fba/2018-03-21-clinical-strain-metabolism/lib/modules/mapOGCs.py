#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annasintsova
"""
import argparse
import os



#takes in ouptut from get-homologs:
# Cluster list:

# ####cluster_list
#format of cluster list:
#
#cluster 5015_thrA size=2 taxa=2 file: 5015_thrA.faa dnafile: 5015_thrA.fna
#: mg1655.gbf
#: HM1.gbk

def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('cluster_list', help = 'list of all the clusters',
                            type=str)
    parser.add_argument('homologs_dir', help = 'directory with cluster files',
                        type=str)
    parser.add_argument('bnum_ref_file', help = 'needed to get bnums for ASAP ids',
                        type=str)
    parser.add_argument('output_file', help = 'name of output csv file',
                        type=str)
    return parser

def getClusterNames(cluster_list):
    ogc = {}
    with open(cluster_list, "r") as fh:
        for line in fh:
            if line.startswith('cluster'):
                name = line.rstrip().split()[1]
                filename = line.rstrip().split()[5]
                ogc[name] = filename
    return ogc

def parseOGCs(cluster_list, homologs_dir):
    #getting all the cluster names and ASAP and PROKKA ids associated with them
    ogc = getClusterNames(cluster_list)
    mapped_ids = {}
    for key, val in ogc.items():
        og = key
        filename = homologs_dir + '/' + val
        asap_id = []
        prokka_id = []# can contain multiple ids
        if os.path.isfile(filename):
            with open(filename, "r") as fh:
                asap = 'NA'
                prokka = 'NA'
                for line in fh:
                    if 'ASAP:' in line:
                        asap = line.split('ASAP:')[1].split(',')[0]
                        asap_id.append(asap)
                    elif 'PROKKA_' in line:
                        prokka = line.split(':')[1].split('|')[0].rstrip()
                        prokka_id.append(prokka)
                mapped_ids[og] = (asap_id, prokka_id)
        else:
            print(filename, ' not found')
            break
    return mapped_ids
# want this intermidiate file as well

def mapToBnum(mapped_ids, bnum_ref_file):
    # mapped id has the following format: cluster_id: ([asap], [prokka])
    id_collection = {}
    for key, value in mapped_ids.items():
        bnum = []
        for asap in value[0]:
            with open(bnum_ref_file) as ref:
                for line in ref:
                    if asap in line:
                        bnum.append(line.split(",")[2])
        id_collection[key] = (value[0], value[1], bnum)
    return id_collection

def output(id_collection, output):
    out_file = open(output, "w+")
    out_file.write("cluster_id\tASAP\tPROKKA\tbnum\n") #tab sepearted
    for key, val in id_collection.items():
        out_file.write("{}\t{}\t{}\t{}\n".format(key, ','.join(val[0]), ','.join(val[1]), ','.join(val[2])))
    out_file.close()


def mapOGCs(cL, homDir, bnum_file, out):
#parse out
    print("Parsing ASAP and PROKKA ids")
    mapped_ids = parseOGCs(cL, homDir)
    print ("Done... Next")
#map to bnums
    print ("Looking for corresponding bnums")
    id_collection = mapToBnum(mapped_ids, bnum_file)
    print ("Done... Next")
#write it output
    print ("Writting out final results")
    output(id_collection, out)
    return "Done"


if __name__ == "__main__":
    args = parser().parse_args()
    cL = args.cluster_list
    homDir = os.path.abspath(args.homologs_dir)
    bnum_file = os.path.abspath(args.bnum_ref_file)
    out = os.path.abspath(args.output_file)
    mapOGCs(cL, homDir, bnum_file, out)


# right now bnum ref is in data/ref/
