# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 09:35:44 2017

@author: yx4vf

#demands,   spectrum,   regenerator,    time (sec)
5:          92,         0,              22
10:         113.5,      0,              127
15:         145.5,      0,              810
"""

from GNmilp1 import *

network_cost = pd.read_csv('networkDT.csv', header=None)/100
network_cost = network_cost.as_matrix()
sn = Network(network_cost)
demands = sn.read_demands('demands_30.csv')
demands = demands.iloc[:15]
demands.source = demands.source-1
demands.destination = demands.destination-1

iteration_history = sn.iterate(demands, mipstart=True, mipfocus=1, timelimit=120, method=2, 
                               mipgap=0.001, outputflag=1, 
                               FeasibilityTol=1e-9, IntFeasTol=1e-9, 
                               OptimalityTol=1e-9)

#model = sn.solve_all(demands)