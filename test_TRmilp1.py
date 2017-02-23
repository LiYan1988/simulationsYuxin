# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:18:35 2017

@author: yx4vf

save solutions from TR and input to GN
"""

from TRmilp1 import *
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
demands = demands.iloc[:15]

iteration_history = sn.iterate(demands, mipstart=True, mipfocus=1, 
                               timelimit=600, method=2, 
                               mipgap=0.001, outputflag=1, 
                               FeasibilityTol=1e-9, IntFeasTol=1e-9, 
                               OptimalityTol=1e-9)

#model = sn.solve_all(demands)