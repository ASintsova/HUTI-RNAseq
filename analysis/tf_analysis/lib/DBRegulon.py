#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 13:54:32 2017

@author: root

Learning how to access RegulonDB via python

source: https://zulko.wordpress.com/2012/08/09/e-colis-regulation-network-from-regulondb-to-python/
"""


path = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/tf_analysis/"

#/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/tf_analysis/
# Create an empty graph
# Create path to the RegulonDB files:
    
    
#path
import networkx as nx
import re

from pylab import *



eColiNetwork = nx.DiGraph(name = 'E.coli Network')

def getRegulonNetwork():
    tfactor_list = []
    tf_reg ={}
    tf_reg_n = {}
    tf_file = open(path + "t_factor_d_tmp.txt")

    for line in tf_file:
        if not line.startswith("#"):
            tf_id = line.rstrip().split("\t")[1]
            tfactor_list.append(tf_id)
       

    reg_file = open(path + "generegulation_tmp.txt")
    for line in reg_file:
        for tf in tfactor_list:
            if tf in line:
                tf_id = line.rstrip().split("\t")[0]
                tf_reg[tf_id] = []


    reg_file = open(path + "generegulation_tmp.txt")
    
    for line in reg_file:
        if not line.startswith("#"):
            tf_u = line.rstrip().split("\t")[0]
            tf_d = line.rstrip().split("\t")[6]
            act = line.rstrip().split("\t")[8]
            tf_u_name = line.rstrip().split("\t")[1]
            tf_d_name = line.rstrip().split("\t")[7]
            
        
            if tf_u_name not in tf_reg_n.keys():
                tf_reg_n[tf_u_name] = [tf_d_name, act]
            else:
                tf_reg_n[tf_u_name] += tf_d_name,act
    #write out regulon file
    regulon_f = open(path+"network_tf_gene.txt", "w+")
    for key in tf_reg_n.keys():
        for i in range(0, int(len(tf_reg_n[key])/2),2):
            if tf_reg_n[key][i+1] == 'activator':
            
                print (key+  "\t" + tf_reg_n[key][i] + "\t+")
                regulon_f.write(key+  "\t" + tf_reg_n[key][i] + "\t+\n")
            elif tf_reg_n[key][i+1] == 'repressor':
                print (key+  "\t" + tf_reg_n[key][i] + "\t-")
                regulon_f.write(key+  "\t" + tf_reg_n[key][i] + "\t-\n")
            else:
                pass



#regulonFile = open(path + "network_tf_gene.txt")

#for line in regulonFile:
#
#    if not line.startswith(('\n', '\t', '#')):
#        g1, g2, sign = line.split('\t')[:3]
#        g1, g2 = g1.lower(), g2.lower()
#        eColiNetwork.add_edge(g1, g2)
#        eColiNetwork[g1][g2]['sign'] = sign
#
#undirected_net = eColiNetwork.to_undirected()
#
#connected_components = nx.connected_components(undirected_net)
#
# 
#print ('Length of the connected components : ', [ len(c) for c in connected_components])
#
#
#hns_regulators = eColiNetwork.successors('hns') + ['hns']
#print ("Genes regulating acs : ", " ".join( acs_regulators))
#
#
#
#component = nx.subgraph(eColiNetwork, hns_regulators)
#from pylab import *
# 
#figure()
#nx.draw(component, with_labels = True)
#show()
#
#


#Trying to get all the transcription factor bnums, will then pull up what can be
#found in DE genes

dat = path + "data/"
fig = path + "figures/"


def TF_bnums(tf_file, bnum_file):
    out = open(dat + "network_tf_bnum.txt", "w+")
    tf = open (dat + tf_file)
    tfs = []
    for line in tf:
        pred = line.rstrip().split("\t")[0]
        suc = line.rstrip().split("\t")[1]
        sign = line.rstrip().split("\t")[2]
        
        bnum = open(dat+bnum_file)
        pred_b = pred
        suc_b = suc
        for st in bnum:
            gene = st.rstrip().split("\t")[2]
            syn = st.rstrip().split("\t")[2]
            if pred in gene or pred in syn:
                pred_b = st.rstrip().split("\t")[7]
            if suc in gene or suc in syn:
                suc_b = st.rstrip().split("\t")[7]
                
            if re.search(r'^b[0-9]', pred_b) != None  and re.search(r'^b[0-9]', suc_b) != None:
                
                #print(pred_b, "\t", suc_b, "\t", sign)
                out.write(pred_b + "\t" + suc_b + "\t" + sign+  "\n")
                if pred_b not in tfs:
                    tfs.append(pred_b)
                
                break
    return tfs        
#tfs = TF_bnums("network_tf_gene.txt", "EcoData.txt")
# have a problem in the original list...

def upregulatedTFs(tfs, de_genes):
    for tf in tfs:
        de_g = open(dat + de_genes)
        for line in de_g:
            
            if line.rstrip().split(",")[3].strip('"') == tf:
                print (line)
        
#upregulatedTFs(tfs, "all_de_genes.csv")
    
#there are clearly some mistakes here, and I will have to go back and manually reannotate this
# but lets see what I can do with this    
    
    
#test_tf = ["b2805", "b2916", "b0076", "b3828", "b3753", "b3512", "b4187", "b2127",
#           "b0069", "b3515", "b0034", "b0487", "b1040"]


test_tf = ["crp", "cspE", "cspD",  "fis",  "leuO", "fucR", "argP", 
           "metR", "xylR", "gadE", "aidB", "mlrA",  "gadX"]
regulonFile = open(dat + "network_tf_gene.txt")

for line in regulonFile:

    if not line.startswith(('\n', '\t', '#')):
        g1, g2, sign = line.split('\t')[:3]
       # g1, g2 = g1.lower(), g2.lower()
        eColiNetwork.add_edge(g1, g2)
        eColiNetwork[g1][g2]['sign'] = sign

undirected_net = eColiNetwork.to_undirected()
#
#connected_components = nx.connected_components(undirected_net)
#
# 
#print ('Length of the connected components : ', [ len(c) for c in connected_components])
#
#


###

#Testing graph drawing


####
test_tf = ["gadE"]

down_regulators = []


for tf in test_tf:
    down_regulators.append(tf) 
    down_regulators += eColiNetwork.successors(tf)
    #down_regulators += eColiNetwork.predecessors(tf)
    component = nx.subgraph(eColiNetwork, down_regulators)
    col = []
    pos=nx.spring_layout(component)
    for e in component.edges():
        col.append('b' if eColiNetwork[e[0]][e[1]]['sign'] == '-\n'else 'r')
       
 
    figure()
    nx.draw(component,pos, edge_color = col, arrows = False, 
            with_labels = True)
    show()



    
###
#Which genes that are regulated by these TFs are also DE 
test_tf = ["crp", "cspE", "cspD",  "fis",  "leuO", "fucR", "argP", 
           "metR", "xylR", "gadE", "aidB", "mlrA",  "gadX"]
down_regulators = []


for tf in test_tf:
    down_regulators.append(tf) 
    down_regulators += eColiNetwork.successors(tf)
    #down_regulators += eColiNetwork.predecessors(tf)
    component = nx.subgraph(eColiNetwork, down_regulators)

re_genes = component.nodes()


def upregulatedTFs2(tfs, de_genes):
    des = []
    for tf in tfs:
  
        de_g = open(dat + de_genes)
        for line in de_g:
            if line.rstrip().split(",")[2].strip('"') == tf:
                print (line)
                des.append((line.rstrip().split(",")[2].strip('"'), float(line.rstrip().split(",")[1])))
    return (des)
x  = upregulatedTFs2(re_genes, "all_de_genes.csv")
x_names = [a[0] for a in x]
x_val =[a[1] for a in x]



down_regulators = ['gadE']

down_regulators += eColiNetwork.successors('gadE')

inc = set(down_regulators).intersection(x_names)

    #down_regulators += eColiNetwork.predecessors(tf)
gad = nx.subgraph(eColiNetwork, down_regulators)
col = []
ncol = []

for e in gad.edges():
    col.append('b' if eColiNetwork[e[0]][e[1]]['sign'] == '-\n'else 'r')

for n in gad.nodes():
    ncol.append('g' if n in inc else 'r')
    
figure()
nx.draw(gad, edge_color = col, node_color = ncol, arrows = False, 
            with_labels = True)

show()


