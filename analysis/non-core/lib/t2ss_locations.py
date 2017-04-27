#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 14:14:32 2017

@author: asintsova
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ori_locs = {'HM1':3631354 ,'HM3':2181548, 'HM6':2822711, 'HM7':3306317, 
            'HM14':3831641, 'HM17':2557959,'HM43':2345263, 
            'HM54': 788250, 'HM56':761996, 'HM57':2827113,
            'HM66':2648967, 'HM68':3298219, 'HM86':686494}

termi_locs = {'HM1':1023758, 'HM3':22243, 'HM6':208333, 'HM7': 807629, 
              'HM14':1425419, 'HM17':5084961,'HM43':4882822, 
              'HM54':3298947, 'HM56':3185894, 'HM57':180836,
            'HM66':131076 , 'HM68':833086, 'HM86':3252017}



all_locs = pd.read_table("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/non-core/data/secretion_ed.gff",
                         header = None)

all_locs.columns = ['genome', 'location', 'gene_name']
x = list(all_locs['gene_name'])
y = [i.split(';')[0].split('=')[1] for i in x] 

locs = all_locs.replace(x, y)
plt.figure(1)

i = 1
for key in ori_locs.keys():
    print(key)

    genome = key
    gm_locs = locs.loc[locs['genome'].str.contains(genome+'_')]
    plt.subplot(7,2,i)
    plt.hist(gm_locs['location'], 50)
    plt.axvline(x=ori_locs[key], color='r', linestyle='dashed', linewidth=2)
    plt.axvline(x=termi_locs[key], color='b', linestyle='dashed', linewidth=2)
    plt.title(genome, fontsize=7)
    plt.ylim(0, 12)
    i+=1
  
#ax0, ax1, ax2, ax3 = axes.flatten()