total_reads_alinged = {'UR86':14606344, 'UTI86':6413793}

import re

def RPKM (counts, gene_info):
	"""
	RPKM = (10^9 * C)/(N * L), with 

C = Number of reads mapped to a gene
N = Total mapped reads in the experiment
L = gene length in base-pairs for a gene

"""
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

	print ('gene\tUrine\tUTI')
	for i in range(len(L)):
	
		print (str(names[i]) + "\t" + str(ur_RPKM[i]) +"\t"+ str(uti_RPKM[i]))
		
x, y = RPKM ("HM86_X_counts", "HM86_papX.df")

print (x)

print (y)
