# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 22:51:13 2017

@author: misWin
"""

import pandas as pd
import numpy as np
nsf = []
idx = 0
with open('nsf-24-map.txt', 'r') as f:
    for line in f:
        if idx==0:
            idx += 1
            continue
        line = line.lstrip('<')
        line = line.rstrip('>,\n')
        line = line.split(',')
        if len(line)==1:
            line = line[0].split()
        line = [int(i) for i in line]
        nsf.append(line)
        idx += 1
        print(line)
        
nsf = pd.DataFrame(nsf)
nsf.columns = ['id', 'source', 'destination', 'length']
nsf.id = nsf.id-1
nsf.source = nsf.source-1
nsf.destination = nsf.destination-1
nsf.length = np.ceil(nsf.length/100)

cost_matrix = np.ones((24, 24))*np.inf
for id in nsf.id:
    cost_matrix[nsf.source[id], nsf.destination[id]] = nsf.length[id]
    
cost_matrix = pd.DataFrame(cost_matrix)
    
cost_matrix.to_csv('nsf-24nodes.csv', header=False, index=False)

cost_matrix = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
