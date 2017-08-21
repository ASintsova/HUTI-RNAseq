#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 14:00:50 2017

@author: root
"""
import argparse
import os


#takes in ouptut from get-homologs: 
# Cluster list:
    
# ####.cluster_list
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
    fh = open(cluster_list)
    for line in fh:
        if line.startswith('cluster'):
            name = line.rstrip().split()[1]
            filename = line.rstrip().split()[5]
            ogc[name] = filename
               
    return ogc

def parseOGCs(cluster_list, homologs_dir, bnum_ref_file, output):
    
    #getting all the cluster names
    ogc = getClusterNames(cluster_list)

    
    
    mapped_ids = {}
    for key, val in ogc.items():
        og = key
        filename = homologs_dir + '/' + val
        #print (filename)
        asap_id = []
        prokka_id = []# can contain multiple ids
        
                    
        if os.path.isfile(filename):
           # print(filename)
            fh = open(filename)
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
    #print(mapped_ids)
    bnum_ids = mapToBnum(mapped_ids, bnum_ref_file)
    out_file = open(output, "w+")
    out_file.write("cluster_id\tmg1655\thm01\n") #tab sepearted
    for key, val in bnum_ids.items():
        out_file.write(key + '\t' +','.join(val[0]) + '\t' + ','.join(val[1]) + '\n')
            
            
def mapToBnum(mapped_ids, bnum_ref_file): # want to replace ASAP ids with bnums

    for key, val in mapped_ids.items():
        asap_to_map = 'NA'
        bnums = []
        prokkas = val[1]
        fh = open(bnum_ref_file) 
        for line in fh:
        
            if line.startswith("ABE"):
                asap_to_map = line.rstrip().split(',')[0]
                bnum = line.rstrip().split(',')[2]

        
            if asap_to_map in val[0]:
                bnums.append(bnum)
                    
        new_val = (bnums, prokkas)
        mapped_ids[key] = new_val     
    return mapped_ids
    


args = parser().parse_args()

cL = args.cluster_list
homDir = os.path.abspath(args.homologs_dir)
out = os.path.abspath(args.output_file)
#print(getClusterNames(args.cluster_list))
bnum_file = os.path.abspath(args.bnum_ref_file)



parseOGCs(cL, homDir, bnum_file, out)
#mapToBnum(bnum_file)








            