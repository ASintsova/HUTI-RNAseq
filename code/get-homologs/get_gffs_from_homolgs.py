# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:50:55 2017

@author: annasintsova
"""
import sys
filename_in = str(sys.argv[1])
filename_out = str(sys.argv[2])

def getGFFs (filename_in, filename_out):
 
    gff = open(filename_out, "w+")
    with open (filename_in, 'r') as fh:
        for line in fh:
            if line.startswith(">>"):
                gene_id = line.rstrip()[2:]
                
                
            elif line.startswith(">"):
                next_line = line
                description = next_line.split("|")
                genome = description[2] + "_1"
                gene_name = description[3]
                start = description[5].split(":")[1].split("-")[0]
                end = description[5].split(":")[1].split("-")[1]
                strand = description[5].split(":")[2].split("^^")[0].strip()
                prokka = description[0].split(":")[1]
                if strand == '1':
                    strand = '+'
                elif strand == '-1':
                    strand = '-'
                else:
                    strand = '.'

                gff.write('{}\tcustom\tCDS\t{}\t{}\t.\t{}\t.\tgene_id "{}"; gene_name "{}"; PROKKA_ID "{}"\n'.format
                      (genome, start, end, strand, gene_id, gene_name, prokka))     
                
getGFFs(filename_in, filename_out)