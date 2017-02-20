# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:11:11 2017

@author: annasintsova
"""
#from "/Users/annasintsova/Desktop/Bioinformatics/Rosalind/findingOrigin" import findingOrigin as fo


#Open genome
import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

path = "/Users/annasintsova/Downloads/RNA_reference_genomes/"


genome_info = {'HM1':(5284090,3631354),#length, origin
'HM3':(4729556,2181548),
'HM6':(5344257,2822711),
'HM7':(4886521,3306317),
'HM14':(5037986,3831641),
'HM17':(5087802,2557959),
'HM43':(4901204,2345263),
'HM54':(5546219,788250),
'HM56':(5040698,761996),
'HM57':(5152486,2827113),
'HM66':(5168961,2648967),
'HM68':(5054821,3298219),
'HM86':(5220431,686494)}

def openMultiFasta(filename):
	
	fhand = open (filename)
	id = ''
	seq = ''
	seq_list = []
	id_list = []
	#extract all sequences and ids and put them together in a dictionary
	for line in fhand:	
		if line.startswith('>'):
			if len(seq) != 0:
				seq_list.append(seq)
			id = line.rstrip()[1:]
			
			seq = ''
			id_list.append(id)
			pass
		else:
			seq += line.rstrip()		
	seq_list.append(seq)
	seq_dict = {}
	for i in range (len(seq_list)):
		seq_dict[id_list[i]]= seq_list[i]
	return seq_dict

geneLoc = 163441
genome = "HM7"
def getGeneLocs(seq_dict, geneLoc, genome):
    origin = genome_info[genome][1]
    print(origin)
    length = genome_info[genome][0]
    print(length)
    contigs = natural_sort([id for id in seq_dict.keys()])
    print (contigs)
    genome_length = 0
    originDist = None
    for seq in contigs:
        print (len(seq_dict[seq]))
        genome_length += len(seq_dict[seq])
        
        if geneLoc < genome_length:
            if geneLoc < origin:
                originDist = 1+((geneLoc - origin)/length)
            else:
                originDist = (geneLoc - origin)/length
        else:
            pass
            
    return originDist
    
    
    