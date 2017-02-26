# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:35:22 2017

@author: yx4vf
"""

from GNmilp1 import *

network_cost = np.array([[INF, INF, INF, INF], 
                         [70, INF, INF, INF],
                         [INF, 70, INF, INF],
                         [70, INF, INF, INF]]).T
sn = Network(network_cost)
demands = pd.DataFrame({'id':[0, 1, 2], 'source':[0, 0, 1], 
                        'destination':[2, 2, 2], 'data_rates':[50, 100, 100],
                        'TR':[75.377, 74.154, 74.154]})
demands = demands[['id', 'source', 'destination', 'data_rates', 'TR']]

iteration_history = sn.iterate(demands, mipfocus=1, timelimit=240, method=2, 
                               mipgap=0.001, outputflag=1, 
                               FeasibilityTol=1e-9, IntFeasTol=1e-9, 
                               OptimalityTol=1e-9)

s = {}
for k in iteration_history[0]['solutions'].keys():
    s[k] = sn.extract_history(iteration_history, k)
    
print(iteration_history[0]['solutions'].keys())