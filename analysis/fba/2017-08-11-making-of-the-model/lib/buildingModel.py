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
    
    for reaction in baseModel.reactions[:]: #go through all reactions in the model, and modify gene_reaction_rules
        print (reaction.id)
        old_rule = reaction.gene_reaction_rule
        old_RR = reaction.gene_reaction_rule.split() #list
        old_RR = [x for x in old_RR if x!= '(' and x != ')' ] # remove brackets
        new_RR = ''
        print(old_RR)
        
        if 'NONE' in old_rule or len(old_RR) < 1:
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
        print('new ', new_RR) 
               
        if 'NA' in new_RR:
            print('deleting')
            reaction.delete()# what are the consequences of deleting reactions, am I creating gaps/problems
                          #in the model?
        else: 
            print('converting')
            reaction.gene_reaction_rule = new_RR
    
    return baseModel
    
#args = parser().parse_args()
##map_file = os.path.abspath(args.map_filename)
#geneMap = getGeneMap(os.path.abspath(args.map_filename))
#
#model = updateReactions(os.path.abspath(args.model_filename), geneMap)
#
#cobra.io.write_sbml_model(model, os.path.abspath(args.new_model_filename))
#





### reading in a sbml file###

#model_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2017-08-11-making-of-the-model/data/mg1655_model.XML"


#bM = cobra.io.read_sbml_model(model_file)

