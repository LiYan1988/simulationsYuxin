# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 11:38:52 2017

@author: yx4vf
"""

#from gurobipy import *
#import numpy as np
#
#np.random.seed(0)
#
#def create_model():
#    model = Model()
#    idx = np.arange(100)
#    tmp = np.random.randint(0, 10, idx.shape)
#    obj = {idx[i]:tmp[i] for i in range(100)}
#    x = model.addVars(idx, vtype=GRB.BINARY, name='x')
#    y = model.addVar(vtype=GRB.BINARY, name='y')
#    model.addSOS(GRB.SOS_TYPE2, x)
#    model.setObjective(x.prod(obj)+y, sense=GRB.MAXIMIZE)
#    
#    return model
#    
#if __name__=='__main__':
#    model = create_model()
#    model.optimize()
#    
##    for x in model.getVars():
##        print(x)
#        
##    y = model.getVarByName('y')
##    y.VType = GRB.INTEGER
##    y.ub = 10
##    y.lb = 0
#    
#    for x in model.getVars():
#        if x.VarName.split('[')[0]=='x':
#            x.ub = 2
#            x.vtype = GRB.INTEGER
#
#    model.optimize()

from milp4 import *

#np.random.seed(0)

#network_cost = pd.read_csv('networkDT.csv', header=None)/100
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost)
n_demands = 50
demands = sn.create_demands(n_demands, modulation='bpsk', low=40, high=100)
demands.to_csv('nsf24-demand1.csv')

idx_iter = 9
idx_stage = 9
n_demands_initial = 5
n_demands_per_stage = 5
demands_added, demands_fixed = sn.scheduler(demands, idx_iter, idx_stage, 
                                            n_demands_initial, 
                                            max_added=n_demands_per_stage)
print(demands_added)
print(demands_fixed)