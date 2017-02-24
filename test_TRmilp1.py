# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:18:35 2017

@author: yx4vf

save solutions from TR and input to GN
#demands,   TR time,    GN time
10,         40.7,       129.8

"""

import TRmilp1 as tr
import GNmilp1 as gn
import numpy as np
import pandas as pd

np.random.seed(0)

network_cost = pd.read_csv('networkDT.csv', header=None)/100
network_cost = network_cost.as_matrix()
sntr = tr.Network(network_cost)
n_demands = 50
demands = sntr.create_demands(n_demands, modulation='bpsk', low=40, high=100)
demands = demands.iloc[:15]

iteration_history_tr = sntr.iterate(demands, mipstart=True, mipfocus=1, 
                                    timelimit=120, method=-1, mipgap=0.001) 

sngn = gn.Network(network_cost)

iteration_history_gn = sngn.iterate(demands, mipstart=True, mipfocus=1, 
                                    timelimit=120, method=-1, mipgap=0.001)