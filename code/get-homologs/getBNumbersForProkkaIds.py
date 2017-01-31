# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 12:12:50 2017

@author: annasintsova

"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:50:55 2017

@author: annasintsova
"""
import sys
import re
#filename_in = str(sys.argv[1])
#filename_out = str(sys.argv[2])
#funct = str(sys.argv[3])
import sys

EC_db_file = "/Users/annasintsova/Documents/TEMP_DATA/EcoData013017-115240.txt"

#for each strain go through old gff file and get prokka_ids, contigs, gene location,
#and UniProt number if there is any
strain_nums = [1, 3, 6, 7, 14, 17, 43, 54, 56, 57, 66, 68, 86]

filenames = []

for strain in strain_nums:
    name = "/Users/annasintsova/Documents/TEMP_DATA/get-homologs-gffs/HM" + str(strain) + "_gh_final.gff"
    filenames.append(name)

    
    
#fnl = "/Users/annasintsova/Documents/TEMP_DATA/RNA_reference_genomes_gff/HM1_edited_gh_FNL.gff"


def SearchForBNumbers(FNL_gff, EC_db_file):
    out_file = "/Users/annasintsova/Documents/TEMP_DATA/get-homologs-gffs/HM" + re.findall("\d+", FNL_gff)[0]+"_bnums.txt"
    fn = open(out_file, "w+")
    fn.write("GENE_ID\tB_NUM")
  
    attributes = []
    fh = open(FNL_gff)
    fh.readline()
    bnums= []
    
    for line in fh:
        
        line = line.rstrip()
        info = line.split("\t")[8].split(";")
        gene_id = info[0].split("=")[1]
        prokka = line.split("\t")[0].split("_")[0]+":"+info[2].split("=")[1]
        uniprot = info[3].split("=")[1]
        product = info[4].split("=")[1]
        symbol = info[1].split("=")[-1].split("_")[0]
        attributes = [symbol, uniprot, product]
        
        bnums = SearchDB(EC_db_file, attributes)
        print (gene_id + "\t" + bnums[0])
        fn.write(gene_id + "\t" + bnums[0] + "\n")
        


def SearchDB(EC_db_file, list_of_attributes):
    
    b_nums = "NA"
    db_descriptor = "NA"
    sanity_check = "NA"
    
    for item in list_of_attributes:
        
        db =open(EC_db_file, 'rt', encoding='latin1')
        for line in db:
            line = line.strip()
            if line.find(item) != -1:
                b_nums = line.split('\t')[7] 
                db_descriptor = line.split('\t')[4]
                sanity_check = line.split('\t')[2]
                break
        if b_nums != "NA":
            break
        
    return [b_nums, db_descriptor, sanity_check]

#for name in filenames:
    #SearchForBNumbers(name, EC_db_file)
    
    
    
    
fns = []

for i in strain_nums:
    n = ("/Users/annasintsova/Documents/TEMP_DATA/get-homologs-gffs/HM" + str(i) + 
        "_bnums.txt")
    fns.append(n)
    
print (fns)    
def BNumQC(filenames):
    fn = "/Users/annasintsova/Documents/TEMP_DATA/get-homologs-gffs/gene_id_to_bnum.txt"
    out_file = open(fn, "w+")
    out_file.write("GENE_ID\tB_NUM")
    id_to_bnum = {}
    for name in filenames:
        print (name)
        fh = open(name)
        fh.readline()
        for line in fh:
            
            ids = line.rstrip().split("\t")[0]
            nums = line.rstrip().split("\t")[1]
           
            if ids not in id_to_bnum:
                id_to_bnum[ids] = [nums]
                #print (id_to_bnum[ids])
            else:
                id_to_bnum[ids] += [nums]
                #print (id_to_bnum[ids])
    for keys in id_to_bnum.keys():
        id_to_bnum[keys] = set(id_to_bnum[keys])
        s = ''
        x = id_to_bnum[keys]        
        for i in x:
            s += str(i) + ";"
        #print (keys + "\t" + s.rstrip(";"))
        out_file.write(keys + "\t" + s.rstrip(";") + "\n")
    return id_to_bnum
        
x = BNumQC(fns)        
        
        