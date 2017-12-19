#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 16:42:01 2017

@author: root
"""


from __future__ import print_function
import cobra
import cobra.test
import os
import argparse
import libsbml

from cobra.util.solver import linear_reaction_coefficients

#tester file

### reading in a sbml file###

model_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/hm01/2017-09-01-hm01-model-2.xml"
model = cobra.io.read_sbml_model(model_file)


model_file2 = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/hm01/hm01_model2.xml"
bM2 = cobra.io.read_sbml_model(model_file2)

model = cobra.test.create_test_model("textbook")

    #update all genes
for g in model.genes:
    #print(g.id)
    if g.id.startswith('b'):
        print(g.name)
        #baseModel.genes.remove(g)
        
        
        
model_file0 = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/mg1655_2011_model.xml"

bM3 = cobra.io.read_sbml_model(model_file0)
print(len(bM3.reactions))

for reaction in bM3.reactions[:50]:
    with bM3 as bM3:
        reaction.knock_out()
        growth_rate = bM3.optimize().objective_value 
                #print('%s blocked (bounds: %s), new growth rate %f' % (reaction.id, str(reaction.bounds), baseModel.objective.value))
                
    if growth_rate > 0:
        reaction.delete()
print(len(bM3.reactions))
      
