#THIS IS NOT WORKING


from __future__ import print_function
import cobra
from cobra.flux_analysis import gapfill

import os
import argparse
import re

#read in the model
print("getting base model")
base = cobra.io.read_sbml_model("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data/ref/mg1655_2011_model.xml")
print("done")

exchange_R = cobra.Model("exchange reactions")

model = cobra.io.read_sbml_model("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data/HM1_homologues/HM1_test.xml")

for i in [i.id for i in model.reactions if 'EX' in i.id]:
    r = model.reactions.get_by_id(i)
    r.delete()


for i in [i.id for i in base.reactions if 'EX' in i.id]:

    print("Adding reaction to exchange model: {}".format(i))
    reaction = base.reactions.get_by_id(i)
    reaction.lower_bound = -1000.0
    reaction.upper_bound = 0.0
    print(reaction in model.reactions)
    exchange_R.add_reaction(reaction.copy())

ts = exchange_R.reactions[2]
print(ts)

#print(exchange_R.reactions[300].lower_bound)
print(ts in model.reactions)
solution = model.optimize()
if solution.f < 0.00001:
    print("gap filling")
    #result = gapfill(model, base)
    result = gapfill(model, exchange_R, exchange_reactions=True)

    for reaction in result[0]:
        print(reaction.id)

else:
    print("No need for gapfilling")
