#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


author: annasintsova
takes in gff file that has info for all the genomes and list of gene_ids that
you are interested in and returns a file with




1. get specific prokka number for each gene_id

2. get KO

3. get bnum

4. edit


"""


gff_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/new_upec_homologs.gff"

gene_ids_txt = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/data/combined_non_core_de_genes.csv" 

prokka_out = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/data/non_core_with_prokkas.txt" 


def getProkkaForGeneID(gff_file, gene_ids_txt, prokka_out):

    gi = open(gene_ids_txt)
    header = gi.readline()
    gene_ids = []
    for line in gi:
        
        gene_ids.append(line.rstrip().split(",")[0].strip('"'))
    out = open (prokka_out, "w+")

    title = 'gene_id\t'
    p_d = {'HM1':'a','HM3':'b',  'HM7':'d','HM14':'d', 'HM17':'e', 'HM43':'f',
                   'HM54':'g', 'HM56':'h', 'HM66':'i','HM57':'j', 'HM68':'k', 'HM86':'l','HM6':'c'}

    for key in sorted(p_d.keys()):
            title += str(key) + '\t'
    out.write(title.rstrip("\t")+ "\n")

    for gene in gene_ids:
        
        s = gene + '\t'
        fh = open(gff_file)
        prokka_dict = {'HM1':'NA','HM3':'NA',  'HM7':'NA','HM14':'NA', 'HM17':'NA', 'HM43':'NA',
                   'HM54':'NA', 'HM56':'NA', 'HM66':'NA','HM57':'NA', 'HM68':'NA', 'HM86':'NA','HM6':'NA'}
        for line in fh:
            
            if gene in line:
                words = line.rstrip().split("\t")
                genome = words[0].split("_")[0]
                prk = words[8].split("PROKKA_ID=")[1]
                prokka_dict[genome] = prk
        for key in sorted(prokka_dict.keys()):
            s += prokka_dict[key] + '\t'
        out.write(s.rstrip("\t") + "\n")

getProkkaForGeneID(gff_file, gene_ids_txt, prokka_out)



###########
"""

takes in a concatenation of all the kegg files Evan, prokka_out, 
and creates file linking each gene_id to KO +other info

"""

kegg_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/kegg_annot_cmb.txt"


KO_out = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/data/non_core_with_KOs.txt" 

def getKO(kegg_file, prokka_out, KO_out):

    p_d = {'HM1':'a','HM3':'b',  'HM7':'d','HM14':'d', 'HM17':'e', 'HM43':'f',
                   'HM54':'g', 'HM56':'h', 'HM66':'i','HM57':'j', 'HM68':'k', 'HM86':'l','HM6':'c'}

    out = open(KO_out, "w+")
    out.write("gene_id\tKO\tsymbol\terror\n")
    
    fh = open(prokka_out)
    fh.readline()
    for line in fh:
        words = line.rstrip().split("\t")
        gene_id = words[0]
        ko = "NA"
        ko2= "NA"
        sym2 = "NA"
        sym = "NA"
        error = 'No'
        
        i = 1
        for key in sorted(p_d.keys()):
            p_d[key]= words[i]
            i+=1
        prk1 = ''
        prk2 = ''
        for key, value in p_d.items():
            if value != "NA" and prk1 == '':
                prk1 = (key, value)
            elif value != "NA" and prk2 == '':
                prk2 = (key, value)
            elif prk1 != '' and prk2 != '':
                break
            else:
                pass
            
        
        kegg_f = open(kegg_file)
        for line in kegg_f:
            
            if prk1[0]+ "~" in line and prk1[1] in line:
                info = line.rstrip().split("\t")
                if info[4].startswith("K"):
                    ko = info[4]
                sym = info[3].split('(')[0]
                #print(ko)
                #print(sym)
            
            if prk2[0]+ "~" in line and prk2[1] in line:
                info = line.rstrip().split("\t")
                if info[4].startswith("K"):
                    ko2 = info[4]
                sym2 = info[3].split('(')[0]
            
            if ko != ko2 and ko != 'NA' and ko2!= 'NA' :
                error = 'Yes'
            
        print (gene_id, "\t", ko, "\t", sym, "\t", error)
        out.write(gene_id + "\t" + ko + "\t"+ sym+ "\t"+ error+"\n")

getKO(kegg_file, prokka_out, KO_out)



###################



kegg_file = "../data/gene_ids_with_KO.txt"
bnum_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/annotations/kegg/eco_genes.txt"
output_txt = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/data/non_core_gene_info.txt"

def editDEgenelist(KO_out,bnum_file, kegg_file, output_txt ):
    out = open (output_txt, "w+")
    fh = open(KO_out)
    header = fh.readline()
    print (header.rstrip() + "\tbnum\n")
    
    out.write(header.rstrip() + "\tbnum\n")
    for line in fh:
        info = line.rstrip().split("\t")
        gene_id= info[0].strip('"')
        gene_name = info[0].strip('"').split("_")[1]
        ko  = info[1]
        sym = info[2]
        bnum = "NA"
        #print (gene_id)
        
                
        bf = open(bnum_file)
        for row in bf:
            if ko!= "NA" and ko in row:
                bnum = row.rstrip().split()[0].split(":")[1]
            elif bnum =="NA" and sym in row:
                bnum = row.rstrip().split()[0].split(":")[1]
            elif bnum == "NA" and gene_name in row:
                bnum = row.rstrip().split()[0].split(":")[1]
        
        print(line.rstrip() + "\t"  + bnum + "\n")
        out.write(line.rstrip() +  "\t" + bnum + "\n")
        
editDEgenelist(KO_out,bnum_file, kegg_file, output_txt )        
        
        
        


