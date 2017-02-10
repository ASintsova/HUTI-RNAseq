#!/usr/bin/python

import sys



total_reads_alinged = {'UR86':14606344, 'UTI86':6413793, 
'UR1':16480318, 'UR3':20927535, 
'UR6':22847371, 'UR7':20980468, 
'UR14': 21533815, 'UR17':19360289, 
'UR43': 18239818, 'UR54':21162541, 
'UR56':17130845, 'UR57':18966744, 
'UR66': 16736066, 'UR68': 15562708, 
'UTI1': 3717040, 'UTI3':8059075, 
'UTI6':1534247, 'UTI7': 683350, 'UTI14':12968218,
'UTI17': 1842583, 'UTI43':2568762, 
'UTI54': 6301997,'UTI56':14935942, 
'UTI57':301746, 'UTI66':79859, 'UTI68':746901}


import re
import sys

"""
	RPKM = (10^9 * C)/(N * L), with 

C = Number of reads mapped to a gene
N = Total mapped reads in the experiment
L = gene length in base-pairs for a gene



"""
# going to do this seperately for each sample: since gene length between homologs
# might vary
# counts filename == key in total_reads_dictionalry
def RPKM  (counts, gh_gff):   
    fh = counts + ".RPKM"
    newfile = open(fh, 'w+')
    newfile.write("GENE_NAME\t" + counts + "\t\n")
    
    L = []
    names = []
    gene_L = {}
    gene = open(gh_gff)
    gene.readline()
	
    for line in gene:
        line = line.rstrip()
        words = line.split("\t")
        L.append(int(words[4]) - int(words[3]))
        names.append(words[8].split(";")[0].split()[1][1:-1]) ### this is not correct
    if len(names) == len(L):
        for i in range (0, len(names)):
            gene_L[names[i]] = L[i]
   
          
    N =  total_reads_alinged[counts]
           
    read_counts = open(counts)
    
    for line in read_counts:
        words = line.rstrip().split("\t")
        if words[0] not in gene_L:
            pass #Needs to give a warning
        else:
            gene_len = gene_L[words[0]]
        
            RPKM = (float(words[1])*10**9)/(N*gene_len)
            #print(words[0] + "\t" + str(RPKM))
            newfile.write(words[0] + "\t" + str(RPKM)+"\n")

counts = str(sys.argv[1])
gh_gff = str(sys.argv[2])	
RPKM (counts, gh_gff)
#
