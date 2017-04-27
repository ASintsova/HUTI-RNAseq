#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:11:24 2017

@author: aasintsova
"""

gene_ids_with_KO = "../data/DE_genes_withKO_and_gene_names.csv"
#gene_ids_with_KO = "../data/test.csv"

bnum_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/eco_genes.txt"
output_txt = "../data/DE_genes_withKO_and_gene_names.csv"


def getBnums(gene_ids_with_KO, bnum_file, output_txt):
    
    fh = open (gene_ids_with_KO)
    header = fh.readline().rstrip()
    out = open (output_txt, "w+")
    out.write(header + "," + '"bnum"\n')
    print(header + "," + '"bnum"')
    for line in fh:
        bnum = "NA"
        ko = line.rstrip().split(",")[1].strip('"')
        if not ko.startswith("K"):
        
            pass
        else:
            bf = open(bnum_file)
            for l in bf:
                if ko in l:
                    bnum = l.rstrip().split()[0].split(":")[1]
            
                
        print(line.rstrip() + ',"' + str(bnum) + '"\n')
        out.write(line.rstrip() + ',"' + str(bnum) + '"\n')
        
        
getBnums(gene_ids_with_KO, bnum_file, output_txt)