#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:44:12 2017

@author: root
"""

#gene_ids_with_KO = "../data/DE_genes_KO_ids.txt"



gene_ids_with_KO = "../data/DE-seq-results-all-core.csv"
kegg_file = "../data/gene_ids_with_KO.txt"
bnum_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/eco_genes.txt"
output_txt = "../data/DE-seq-results-all-core-edited.csv"

def editDEgenelist(gene_ids_with_KO,bnum_file, kegg_file, output_txt ):
    out = open (output_txt, "w+")
    fh = open(gene_ids_with_KO)
    header = fh.readline()
    out.write(header.rstrip() + ',"sym","KO","bnum"\n')
    for line in fh:
        info = line.rstrip().split(",")
        gene_id= info[0].strip('"')
        gene_name = gene_id.split("_")[1]
        KO = "NA"
        sym = "NA"
        bnum = "NA"
        #print (gene_id)
        
        kf = open(kegg_file)
        
        for l in kf:
            
            if gene_id in l:
                info = l.rstrip().split("\t")
                if info[1].startswith("K"):
                    KO = info[1]
                sym = info[3].split("(")[0]
                
        bf = open(bnum_file)
        for row in bf:
            if KO != "NA" and KO in row:
                bnum = row.rstrip().split()[0].split(":")[1]
            elif bnum =="NA" and sym in row:
                bnum = row.rstrip().split()[0].split(":")[1]
            elif bnum == "NA" and gene_name in row:
                bnum = row.rstrip().split()[0].split(":")[1]
        
        print(line.rstrip() + "," + sym + "," + KO + "," + bnum + "\n")
        out.write(line.rstrip() + "," + sym + "," + KO + "," + bnum + "\n")
        
editDEgenelist(gene_ids_with_KO,bnum_file, kegg_file, output_txt )        
        
        
        