from __future__ import print_function
import cobra

import os
import argparse
import re

bnum = re.compile(r'b[0-9]*')
OGC_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data/HM1_homologues/HM1_OGC.tab"
with open(OGC_file, "r") as fh:
    text = fh.read()
    strain_genes = bnum.findall(text)

#df = pandas.read_table(OGC_file)
#strain_genes = list(df['bnum'])
print (len(strain_genes))

base = cobra.io.read_sbml_model("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data/ref/mg1655_2011_model.xml")
model_genes = list(base.genes)
#x = str(model_genes[1])
#print (x)
#print (x in strain_genes)
genes_to_remove = [g for g in model_genes if str(g) not in strain_genes]
print(len(genes_to_remove))

if len(genes_to_remove)>0:
    cobra.manipulation.delete_model_genes(base, genes_to_remove)
else:
    base._trimmed_genes=[]
    base._trimmed_reactions=[]

print("Deleted: {} genes, {} reactions".format(len(base._trimmed_genes), len(base._trimmed_reactions)))
