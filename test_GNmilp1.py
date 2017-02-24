# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 09:35:44 2017

@author: yx4vf

#demands,   time (sec)
5:          22
10:         127
15:         810
20:         4460, timelimit=600
"""

from GNmilp1 import *
np.random.seed(0)

network_cost = pd.read_csv('networkDT.csv', header=None)/100
network_cost = network_cost.as_matrix()
sn = Network(network_cost)
# read from matlab data
#demands = sn.read_demands('demands_30.csv')
#demands = demands.iloc[:15]
#demands.source = demands.source-1
#demands.destination = demands.destination-1
n_demands = 50
demands = sn.create_demands(n_demands, modulation='bpsk', low=40, high=100)
demands = demands.iloc[:10]

iteration_history = sn.iterate(demands, mipstart=True, mipfocus=1, 
                               timelimit=120, method=-1, 
                               mipgap=0.001, outputflag=1, 
                               FeasibilityTol=1e-7, IntFeasTol=1e-7, 
                               OptimalityTol=1e-7)

#model = sn.solve_all(demands)