#!/usr/bin/python

import sys



total_reads_alinged = {'UR86':14606344, 'UTI86':6413793, 'UR1':16480318, 'UR3':20927535, 
'UR6':22847371, 'UR7':20980468, 'UR14': 21533815, 'UR17':19360289, 'UR43': 18239818, 'UR54':21162541, 
'UR56':17130845, 'UR57':18966744, 'UR66': 16736066, 'UR68': 15562708, 'UTI1': 3717040, 'UTI3':8059075, 'UTI6':1534247, 'UTI7': 683350, 'UTI14':12968218,
'UTI17': 1842583, 'UTI43':2568762, 'UTI54': 6301997,'UTI56':14935942, 'UTI57':301746, 'UTI66':79859, 'UTI68':746901}


import re

def RPKM (counts, gene_info):
	"""
	RPKM = (10^9 * C)/(N * L), with 

C = Number of reads mapped to a gene
N = Total mapped reads in the experiment
L = gene length in base-pairs for a gene

"""
	fh = counts + ".RPKM"
	newfile = open(fh, 'w+')
	ur_RPKM = []
	uti_RPKM = []
	L = []
	names = []
	gene = open(gene_info)
	gene.readline()
	
	#for each gene calculate length
	for line in gene:
		line = line.rstrip()
		words = line.split()
		L.append(int(words[2]) - int(words[1]))
		names.append(words[4])
	
	read_counts = open(counts)
	urine = []
	UTI =  []
	read_counts.readline()
	#for each condition retrieve how many reads mapped
	UR_N = total_reads_alinged['UR' + re.findall("[0-9]+", counts)[0]]
	UTI_N = total_reads_alinged['UTI' + re.findall("[0-9]+", counts)[0]]
	#collect of raw reads for each condition in different list
	for line in read_counts:
		line = line.rstrip()
		w = line.split(',')
		urine.append(int(w[1]))
		UTI.append(int(w[2]))
		
	if len(L) == len(urine) and len(L) == len (UTI):
		for i in range(len(L)):
			ur_RPKM.append ((10**9 * urine[i])/(float(UR_N) * L[i]))
			uti_RPKM.append ((10**9 * UTI[i])/(float(UTI_N) * L[i]))
		#return  ur_RPKM, uti_RPKM
		
	else: 
		print ('ohoh')
	newfile.write('gene\tUrine\tUTI\n')
	print ('gene\tUrine\tUTI')
	for i in range(len(L)):
	
		print (str(names[i]) + "\t" + str(ur_RPKM[i]) +"\t"+ str(uti_RPKM[i]))
		newfile.write(str(names[i]) + "\t" + str(ur_RPKM[i]) +"\t"+ str(uti_RPKM[i])+ "\n")
counts = str(sys.argv[1])
gene_info = filename = str(sys.argv[2])	
RPKM (counts, gene_info)

