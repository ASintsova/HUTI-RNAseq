# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 17:19:47 2016

@author: annasintsova
"""
import re
import sys

#gene list with locations to be identified: first col genome,
#second and third locations within genome, fifth is the LS-BSR identifier
#first line header
filename = str(sys.argv[1])

fh = open(filename)

core_genes_locs = {}
core_genes_names = {}
fh.readline()
for line in fh:
    words = line.rstrip().split()
    core_genes_locs[words[4]] = (words[1], words[2])
    core_genes_names[words[4]] = []
    
#annotation file
ann_file = str(sys.argv[2])
anh = open(ann_file)

for gene in core_genes_locs:
    anh = open(ann_file)
    sear = core_genes_locs[gene]
    location = sear[0]+'|'+sear[1]    
    for line in anh:
        m = re.search(location, line)
    	if m != None:
        	words = line.rstrip().split("\t")
        	description = words[8]
        	if "gene=" in description:
        		
        		core_genes_names[gene] += [description.split("gene=")[1].split(";")[0]]
        	else:
        		core_genes_names[gene] += ["NA"]
        	core_genes_names[gene] += [description.split("product=")[-1]]
        	#print (core_genes[gene])
output = str(filename[:-3] + '.id')
identified = open(output, "w+")

for keys, values in core_genes_names.items():
	s = ''
	for val in values:
		s += val + "\t"
	identified.write(str(keys)+ '\t' + s.rstrip()+ "\n")
    
    
       