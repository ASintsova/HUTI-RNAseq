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


































#def getGFFs (filename_in, filename_out):
# 
#    gff = open(filename_out, "w+")
#    with open (filename_in, 'r') as fh:
#        for line in fh:
#            if line.startswith(">>"):
#                gene_id = line.rstrip()[2:]
#                
#                
#            elif line.startswith(">"):
#                next_line = line
#                description = next_line.split("|")
#                genome = description[2] + "_1"
#                gene_name = description[3]
#                start = description[5].split(":")[1].split("-")[0]
#                end = description[5].split(":")[1].split("-")[1]
#                strand = description[5].split(":")[2].split("^^")[0].strip()
#                prokka = description[0].split(":")[1]
#                if strand == '1':
#                    strand = '+'
#                elif strand == '-1':
#                    strand = '-'
#                else:
#                    strand = '.'
#
#                gff.write('{}\tcustom\tCDS\t{}\t{}\t.\t{}\t.\tgene_id={}; gene_name={}; PROKKA_ID={}\n'.format
#                      (genome, start, end, strand, gene_id, gene_name, prokka))     
#                #print ('{}\tcustom\tCDS\t{}\t{}\t.\t{}\t.\tgene_id "{}"; gene_name "{}"; PROKKA_ID "{}"\n'.format
#                      #(genome, start, end, strand, gene_id, gene_name, prokka))
#
#
#def crossRefWithProkka (gh_gff, prokka_gff):
#    prokka = open(prokka_gff)
#    wrong_locs = []
#    wl = gh_gff[:-4] + "_wrong_locs.txt"
#    
#    #wl_fn = open(wl, "w+")
#    
#    gh = open(gh_gff)
#    prokka.readline()# skip first line
#    prokka_ids = {}
#    
#    fn = gh_gff[:-4] + "_final.gff"
#    #final_gff = open(fn, "w+")
#    #final_gff.write("#gff-veresion 3")
#    
#    UniProt = 'NO_UniProt'
#    for line in prokka:
#        if line.startswith("##"):
#            pass
#    
#        else:
#            words = line.rstrip().split("\t")
#            info = words[8].split(';')
#            
#            prokka_id = info[0].split("=")[1]
#            for i in info:
#                
#                if i.find("UniProt") != -1:
#                    UniProt = i.split(":")[-1]
#            
#            prokka_ids_contig[prokka_id]= [words[0], UniProt]
#            
#            prokka_ids_locs[prokka_id] = (int(words[3]), int(words[4]))
#            
#    for line in gh:
#         
#        words = line.rstrip().rsplit("\t")
#        
#        
#        
#        ID = words[8].split("PROKKA_ID=")[1]
#     
#        new_line = "\t".join(words[1:9])
#        
#        #print (prokka_ids_locs[id])
#        #print(id)
#        if ID not in prokka_ids_locs.keys():
#            print (ID + " not found")
#        elif (int(words[3]), int(words[4])) != prokka_ids_locs[id]:
#            wrong_locs.append(ID)
#            print (prokka_ids_locs[ID])
#        print (str(prokka_ids_contig[ID][0], "\t", new_line+";" + 
#        #final_gff.write(str(prokka_ids_contig[id][0])+"\t" +new_line + "\n")
#       
#    #wl_fn.write(str(wrong_locs))
#    return wrong_locs   
#
#
#gh_gff = "/Users/annasintsova/Documents/TEMP_DATA/Annotations/edited-get-homologs-gff/HM43_edited_gh.gff"
#prokka_gff = "/Users/annasintsova/Documents/TEMP_DATA/Annotations/HM43.gff"


#if funct == '1':
#    
#    getGFFs(filename_in, filename_out) 
#    
#elif funct == '2':
#    
#    crossRefWithProkka (filename_in, filename_out)
#else:
#    print ("Invalid Input")
#        