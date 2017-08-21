# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 11:29:24 2016

@author: annasintsova
"""

#!/usr/bin/python

import sys

def blastFileandEdittoDF (filename , newbedfile):

    strand = ''
    h1 = ''
    h2 = ''
    count = 1
    fhand = open (filename) #filename is name of balst output file
    for line in fhand:
    	if len(line) == 0:
    		break
        line = line.rstrip()
        word = line.split("\t")
        if int(word[2]) < int(word[1]):
            strand = '-'
            h1 = word[2]
            h2 = word[1]
        else:
            strand = '+'
            h1 = word[1]
            h2 = word[2]
        print (word[0]+"\t"+h1+"\t"+h2+"\t"+str(strand)+'\t'+ 'gene_'+str(count) + ';'+ word[4]) 
        newbedfile.write(word[0]+"\t"+h1+"\t"+h2+"\t"+str(strand)+'\t'+ 'gene_'+str(count) + ';'+ word[4]+ "\n")
        count += 1
filename = str(sys.argv[1])
fh = filename + ".df"
newbedfile = open(fh, 'w+')
newbedfile.write("seqnames\tstart\tend\tstrand\tfeature\n")
blastFileandEdittoDF(filename, newbedfile)
