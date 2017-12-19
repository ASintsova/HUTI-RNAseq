#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Aug 15 10:08:05 2017

@author: annasintsova
"""

from __future__ import print_function
import cobra
import cobra.test
import os
import argparse
import libsbml
import random

from corda import reaction_confidence

from corda import CORDA

def parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('map_filename', help = 'file with model ids mapped to prokka ',
                            type=str)
    parser.add_argument('model_filename', help = 'file with base model',
                        type=str)
    parser.add_argument('new_model_filename', help = 'output filename',
                        type=str)

    return parser


def getGeneMap(map_filename): #file mapping model genes to genome of interest genes, tab seperated

    geneMap = {}
    fh = open(map_filename)
    fh.readline()
    for line in fh:
        bnum = line.rstrip().split('\t')[1].split(',')
        prokka = line.rstrip().split('\t')[2].split(',')
        
        for num in bnum:
            
            geneMap[num] = prokka # all the prokka ids map to each of the bnums
        
    return geneMap


def updateReactions(baseModel_file, geneMap):
    
    baseModel = cobra.io.read_sbml_model(baseModel_file)# read in base model
    
    #update all the reactions
    
    for reaction in baseModel.reactions: #go through all reactions in the model, and modify gene_reaction_rules
        print (reaction.id)
        old_rule = reaction.gene_reaction_rule
        print(old_rule)
        old_RR = reaction.gene_reaction_rule.split()
        #list
        old_RR = [x.strip('(').strip(')') for x in old_RR if x!= '(' and x != ')' ] # remove brackets
        print(old_RR)
        new_RR = ''
        if 'NONE' in old_rule or len(old_RR) < 1: #either empty or listed as NONE
            #print(old_RR)
            new_RR = old_rule
        
        elif len(old_RR) == 1:
        
            try:
                new_RR = "( " + ' or '.join(geneMap[old_RR[0]]) + " )" # for long logicals need brackets
            except KeyError:
                new_RR = 'NA'
            
        elif len(old_RR) > 1:
        
            for i in range(0,  len(old_RR)):
            
                if i%2 == 0:
                    try:
                        new_RR += "( "
                        new_RR += ' or '.join(geneMap[old_RR[i]])
                        new_RR += " )"
                    except KeyError:
                        new_RR = 'NA'
                    
                else:
                    new_RR += ' '
                    new_RR += old_RR[i]
                    new_RR += ' '
        #print('new ', new_RR) 
               
        if 'NA' in new_RR:
            print(reaction.id)
            with baseModel as baseModel:
                reaction.knock_out()
                growth_rate = baseModel.optimize().objective_value 
                #print('%s blocked (bounds: %s), new growth rate %f' % (reaction.id, str(reaction.bounds), baseModel.objective.value))
                
            if growth_rate > 0:
                print('deleting ', reaction.id)
                reaction.delete()
        else: 
            print('converting ', old_RR, " to ", new_RR)
            reaction.gene_reaction_rule = new_RR
    #update all genes
    gene_ids = [g.id for g in baseModel.genes]
    for g in gene_ids:
        if g.startswith('PROKKA'):
            pass
        else:    
            G = baseModel.genes.get_by_id(g)
            baseModel.genes.remove(G)
            
    
    return baseModel
    
#args = parser().parse_args()
#map_file = os.path.abspath(args.map_filename)
#geneMap = getGeneMap(os.path.abspath(args.map_filename))
#
#model = updateReactions(os.path.abspath(args.model_filename), geneMap)
#
#cobra.io.write_sbml_model(model, os.path.abspath(args.new_model_filename))
#

gene_ids_with_prokkas_tab = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/gene_ids_with_prokkas.txt"
confScores_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/hm01/uti1_confScore.csv"
confScore_out_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/hm01/uti1_confScore_final.csv"




def mapIDs(gene_ids_with_prokkas_tab, confScores_file, confScore_out_file):
    Sc = open(confScores_file)
    new_Sc = open (confScore_out_file, "w+")
    header = Sc.readline().rstrip() + ',prokka_id\n'
    new_Sc.write(header)
    
    for line in Sc:
        
        gene_name = line.split(',')[0]
        
        ids_with_p = open(gene_ids_with_prokkas_tab)
        for l in ids_with_p:
            if gene_name in l:
                prokka_id = l.split('\t')[1]
                
                new_Sc.write(line.rstrip() + ',' + prokka_id + '\n')
                #print(line.rstrip() + ',' + prokka_id + '\n')
                break
                    
#
#model = cobra.test.create_test_model("textbook")
#options = [-1, 0, 1,2,3]
#gene_conf = {}
#for gene in model.genes:
#    gene_conf[gene.id] = random.choice(options) # needs to be gene ids
#
#rx_conf = {}
#
#for reaction in model.reactions:
#    rule = reaction.gene_reaction_rule
#    if rule != 'NONE' or len(rule) > 0:
#            print('mapping')
#            rx_conf[reaction.id] = reaction_confidence(rule, gene_conf)
#           # 
#    
#opt = CORDA(model, rx_conf)
#opt.build()
#print(opt)
#
#
#model_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/hm01/2017-09-01-hm01-model-2.xml"
#model = cobra.io.read_sbml_model(model_file)
#



def mapReactionConfidences(geneConfScores_csv, model):
    
    cS= open(geneConfScores_csv)
    
    header = cS.readline().strip().split(',')
    print(header)
    scoreInx = header.index('confScore')
    idInx = header.index('prokka_id')
    gene_conf = {}
    for line in cS:
        gene_conf[line.strip().split(',')[idInx]] = int(line.strip().split(',')[scoreInx])
    print(gene_conf)    
    rx_conf = {}   
    for reaction in model.reactions:
        rule = reaction.gene_reaction_rule
        if rule != 'NONE' or len(rule) > 0:
            rx_conf[reaction.id] = reaction_confidence(rule, gene_conf)
    
    
    
    opt = CORDA(model, rx_conf)
    opt.build()
    
    #print(opt)
    return opt
    
    
    

### reading in a sbml file###

#model_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/mg1655_model.XML"


#bM = cobra.io.read_sbml_model(model_file)
