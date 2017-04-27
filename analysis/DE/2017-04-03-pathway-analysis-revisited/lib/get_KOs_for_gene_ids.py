#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:43:32 2017

@author: annasintsova

takes in a concatenation of all the kegg files Evan, gene_ids_with_prokkas (output
from get_prokkas_for_gene_ids.py), and creates file linking each gene_id to KO +
other info

"""

kegg_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/kegg_annot_cmb.txt"

gene_ids_txt = "data/gene_ids_with_prokkas.txt"

output_txt = "data/gene_ids_with_KO.txt"

def getKO(kegg_file, gene_ids_txt, output_txt):

    p_d = {'HM1':'a','HM3':'b',  'HM7':'d','HM14':'d', 'HM17':'e', 'HM43':'f',
                   'HM54':'g', 'HM56':'h', 'HM66':'i','HM57':'j', 'HM68':'k', 'HM86':'l','HM6':'c'}

    out = open(output_txt, "w+")
    fh = open(gene_ids_txt)
    fh.readline()
    #for line in fh:
    #line = fh.readline()
    for line in fh:
        words = line.rstrip().split("\t")
        gene_id = words[0]
        i = 1
        for key in sorted(p_d.keys()):
            p_d[key]= words[i]
            i+=1


        kegg_f = open(kegg_file)
        x = 0
        all_KO = 'NA'
        for line in kegg_f:
            test = ''

            for key, value in p_d.items():
                genome = key +"~"

                if x == 0 and genome in line and value in line:
                    info = line.rstrip().split("\t")
                    KO = info[4]
                    code = info[2]
                    description = info[3]
                    x = 1

               #print (KO)

                if x != 0 and genome in line and value in line:

                    info = line.rstrip().split("\t")
                    test = info[4]

                    if test != KO:

                        if all_KO == 'NA':
                            all_KO = KO + ","+test
                        else:
                            all_KO += ','+test
                        if not KO.startswith('K') and test.startswith('K'):
                            KO = test


        #print (gene_id, "\t", KO, "\t", code, "\t", description,"\t", all_KO)
        out.write(gene_id + "\t" + KO + "\t"+ code+ "\t"+ description+"\t"+ all_KO+ "\n")

getKO(kegg_file, gene_ids_txt, output_txt)
