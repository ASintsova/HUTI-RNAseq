# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 14:52:42 2017

@author: annasintsova

editing Log2FC file to extract gene names where possible

"""

fh = open ("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/2017-01-20-core-genes-DE.csv")
edited_file = open("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/get-homologs/counts/2017-01-23-core-genes-DE-edited.csv", "w+")

edited_file.write(fh.readline().rstrip() + "," + "geneName"+ "\n")

for line  in fh:
    
    words = line.rstrip().split(",")
    
    name = words[0].split("_")

    if len(name[1]) == 4:
        ed_name = name[1] 
    else:    
    
        ed_name = "_".join(words[0].split("_")[1:])
    
    #Make sure weird quotes are not attached
    
    if ed_name[-1] == '"':
        
        ed_name = ed_name[:-1]
    
    # Make sure all start with lower letter
    fm_letter = [c for c in ed_name[1:]]
    new_letter = [ed_name[0].lower()]   
    ed_name = ''.join((new_letter + fm_letter))
    
    print (line.rstrip() + "," + ed_name)
    edited_file.write(line.rstrip() + "," + ed_name + "\n")