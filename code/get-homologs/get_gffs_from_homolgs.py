# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:50:55 2017

@author: annasintsova
"""
import sys
filename_in = str(sys.argv[1])
filename_out = str(sys.argv[2])
funct = str(sys.argv[3])

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

                gff.write('{}\tcustom\tCDS\t{}\t{}\t.\t{}\t.\tgene_id={}; gene_name={}; PROKKA_ID={}\n'.format
                      (genome, start, end, strand, gene_id, gene_name, prokka))     
                #print ('{}\tcustom\tCDS\t{}\t{}\t.\t{}\t.\tgene_id "{}"; gene_name "{}"; PROKKA_ID "{}"\n'.format
                      #(genome, start, end, strand, gene_id, gene_name, prokka))


def crossRefWithProkka (gh_gff, prokka_gff):
    prokka = open(prokka_gff)
    wrong_locs = []
    wl = gh_gff[:-4] + "_wrong_locs.txt"
    
    wl_fn = open(wl, "w+")
    
    gh = open(gh_gff)
    prokka.readline()# skip first line
    prokka_ids_contig = {}
    prokka_ids_locs = {}
    fn = gh_gff[:-4] + "_final.gff"
    final_gff = open(fn, "w+")
    final_gff.write("#gff-veresion 3")
    for line in prokka:
        if line.startswith("##"):
            break
    
        else:
            words = line.rstrip().split("\t")
            
            prokka_id = words[8].split(";")[0].split("=")[1]
            prokka_ids_contig[prokka_id]= words[0]
            prokka_ids_locs[prokka_id] = (int(words[3]), int(words[4]))
   
    for line in gh:
         
        words = line.rstrip().rsplit("\t")
        
        
        
        id = words[8].split("PROKKA_ID=")[1]
     
        new_line = "\t".join(words[1:9])
        
       
        if (int(words[3]), int(words[4])) != prokka_ids_locs[id]:
            wrong_locs.append(id)
        
        final_gff.write(str(prokka_ids_contig[id])+"\t" +new_line + "\n")
       
    wl_fn.write(str(wrong_locs))
    
if funct == '1':
    
    getGFFs(filename_in, filename_out) 
    
elif funct == '2':
    
    crossRefWithProkka (filename_in, filename_out)
else:
    print ("Invalid Input")
        