import os
import pandas as pd
deleted_genes = {}
essent_genes = {}
for root, dirs,files in os.walk('./data', topdown =True):
    for fh in files:
       if "deleted_NOT_FOUND" in fh:
           with open(os.path.join(root, fh)) as fi:
               for line in fi:
                   line = line.rstrip("\n")
                   deleted_genes[line] = deleted_genes.get(line, 0) + 1
       elif "retained_NOT_FOUND" in fh:
           with open(os.path.join(root, fh)) as fi:
               for line in fi:
                   line = line.rstrip("\n")
                   essent_genes[line] = essent_genes.get(line, 0) + 1

bnums = pd.Series(deleted_genes,index=deleted_genes.keys())
bnums.sort(ascending=False)
#print(bnums)

essent_bnums = pd.Series(essent_genes, index=essent_genes.keys())
essent_bnums.sort(ascending=False)
print(essent_bnums)

