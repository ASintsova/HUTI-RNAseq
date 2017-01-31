# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:50:55 2017

@author: annasintsova
"""
import sys
prokka_gff = str(sys.argv[1])
gh_gff = str(sys.argv[2])
#funct = str(sys.argv[3])


    

def getProkkaIds(prokka_gff):
    prokka = open(prokka_gff)  
 
    prokka.readline()
    prokka_ids = {}   
     
    UniProt = 'NA'
    product = "NA"
    for line in prokka:
        if line.startswith("HM"):
                
            words = line.rstrip().split("\t")
          
            info = words[8].split(';')
          
            
            prokka_id = info[0].split("=")[1]
            
            product = words[8].split("=")[-1]
            for i in info: 
                if i.find("UniProt") != -1:
                    UniProt = i.split(":")[-1]
         
            prokka_ids[prokka_id]= [words[0], int(words[3]), int(words[4]), UniProt, product]
    
    return prokka_ids  
            



def CleanGhGff(gh_gff, prokka_ids):
    
    wl = gh_gff[:-4] + "_wrong_locs.txt"
    wl_file = open(wl, "w+")
    wrong_locs = []
    
    fn = gh_gff[:-4] + "_final.gff" 
    
    f_out = open(fn, "w+")
    f_out.write("##gff-version 3\n")
    
    gh = open(gh_gff)
    for line in gh:
        line = line.rstrip()
        gene_info = line.split("\t")[8]
        
        # get prokka id, and look up appropriate contig from the dictionary     
        p_id = gene_info.split("PROKKA_ID=")[1]
        new_contig = prokka_ids[p_id][0]
        
        #make sure gene locations match between two files
        old_locations = (int(line.split("\t")[3]), int(line.split("\t")[4]) )   
        new_locations = (prokka_ids[p_id][1], prokka_ids[p_id][2])
        if  old_locations != new_locations:
            wrong_locs.append(p_id)
        print (str(new_contig)+ "\t"+ "\t".join(line.split("\t")[1:8])+
                "\t"+ str(gene_info)+ "; UniProt=" + str(prokka_ids[p_id][3]) +
                "; product=" + str(prokka_ids[p_id][4]))  
                
        f_out.write (str(new_contig)+ "\t"+ "\t".join(line.split("\t")[1:8])+
                "\t"+ str(gene_info)+ "; UniProt=" + str(prokka_ids[p_id][3]) +
                "; product=" + str(prokka_ids[p_id][4])+ "\n")
    wl_file.write(str(wrong_locs))

def editGhGff(prokka_gff, gh_gff):
    
    P = getProkkaIds(prokka_gff)
    CleanGhGff(gh_gff, P)
    
editGhGff(prokka_gff, gh_gff)






























