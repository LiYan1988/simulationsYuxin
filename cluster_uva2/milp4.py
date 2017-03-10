# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 17:35:05 2017

@author: yx4vf
"""
########################### DO NOT CHANGE!!! #################################
import numpy as np
import pandas as pd
from gurobipy import *
import time
import copy
import pickle
import gc
from decimal import Decimal

# fiber parameters
INF = np.inf # infinity
G = 12.5 # guardband
cofase = 23.86 # ASE coefficient
rou = 2.11*10**-3
miu = 1.705

# objective weight
Nmax = 10 # max number of regenerator circuits per regenerator node
epsilon_total = 1
epsilon_nnn = 0/(Nmax*10)

# modelling parameters
bigM1 = 10**4
bigM2 = 10**4
bigM3 = 2*10**4

# scheduler parameters
n_demands_per_stage = 5
n_iter_per_stage = 5 # 10
timelimit_baseline = 150 # 960
timelimit0 = 20 # 60
time_factor = 1.5
##############################################################################

class Network(object):
    '''The network 
    '''
    
    def __init__(self, cost_matrix, modulation='bpsk'):
        '''Initialize the network topology'''
        self.cost_matrix = cost_matrix
        self.connectivity_matrix = 1-np.isinf(self.cost_matrix)
        self.n_nodes = self.cost_matrix.shape[0]
        self.nodes = list(range(self.n_nodes))
        self.node_pairs = [(i, j) for i in range(self.n_nodes) 
            for j in range(self.n_nodes) if i!=j]
        self.n_node_pairs = len(self.node_pairs)
        links = [[i, j, self.cost_matrix[i,j]] 
            for i in range(self.n_nodes)
            for j in range(self.n_nodes)
            if i!=j and 1-np.isinf(self.cost_matrix[i,j])]
        ids = [i for i in range(len(links))]
        src = [links[i][0] for i in range(len(links))]
        dst = [links[i][1] for i in range(len(links))]
        length = [int(links[i][2]) for i in range(len(links))]
        self.n_links = len(links)
        self.links = pd.DataFrame({'id': ids, 'source': src, 
                                   'destination': dst, 'length':length})
        self.links = self.links[['id', 'source', 'destination', 'length']]
        if modulation=='bpsk':
            self.Noise = (10**4)/3.52 # bpsk
        elif modulation=='qpsk':
            self.Noise = (10**4)/7.03 # qpsk
    
    def create_demands(self, n_demands, modulation='bpsk', distribution='uniform', 
                       low=30, high=400):
        '''Create demands
        '''
        demand_pairs_id = np.random.choice(self.n_node_pairs, n_demands)
        demand_pairs_id = sorted(demand_pairs_id)
        demand_pairs = [self.node_pairs[i] for i in demand_pairs_id]
        if distribution=='uniform':
        # uniform distribution
            data_rates = np.random.randint(low, high, n_demands)
        elif distribution=='normal':
        # normal distribution
            data_rates = np.random.normal(loc=(low+high)/2, scale=(low+high)/4,
                                          size=n_demands)
            data_rates = [int(max(low, data_rates[i])) for i in range(n_demands)]
            data_rates = [int(min(high, data_rates[i])) for i in range(n_demands)]
        src = [x[0] for x in demand_pairs]
        dst = [x[1] for x in demand_pairs]
        ids = [i for i in range(n_demands)]
        demands = pd.DataFrame({'id': ids, 'source':src, 'destination':dst, 
                                'data_rates': data_rates})
        demands = demands[['id', 'source', 'destination', 'data_rates']]

        # choose modulation format
        if modulation=='qpsk':
            qpsk_tr = pd.read_csv('qpsk_TR.csv', header=None)
            qpsk_tr.columns = ['data_rate', 'distance']
            qpsk_tr.distance = qpsk_tr.distance/100
            qpsk_tr.set_index('data_rate', inplace=True)
            tr = [float(qpsk_tr.loc[int(np.round(data_rates[i]))]) 
                for i in range(n_demands)]
        elif modulation=='bpsk':
            bpsk_tr = pd.read_csv('bpsk_TR.csv', header=None)
            bpsk_tr.columns = ['data_rate', 'distance']
            bpsk_tr.distance = bpsk_tr.distance/100
            bpsk_tr.set_index('data_rate', inplace=True)
            tr = [float(bpsk_tr.loc[int(np.round(data_rates[i]))]) 
                for i in range(n_demands)]
            
        demands['TR'] = tr

        return demands
    
    def create_model_gn(self, demands, **kwargs):
        '''Create MILP for GN at the first iteration of a stage, no demands are
        fixed, do not solve the model
        '''
        n_demands = demands.shape[0]
        
        # define supply 
        supply = np.zeros((self.n_nodes, n_demands))
        for n in range(self.n_nodes):
            for d in demands.id:
                if demands.source[d]==n:
                    supply[n, d] = -1

        for n in range(self.n_nodes):
            for d in range(n_demands):
                if demands.destination[d]==n:
                    supply[n, d] = 1
        
        model = Model('GN_{}'.format(demands.shape[0]))
        model.Params.UpdateMode = 1
        
        # define variables
        U = {} # U[a, b] = UsageL[a,b]*Ynode[a,b]
        for l in self.links.id:
            for d in demands.id:
                U[l, d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                 name='U_{}_{}'.format(l, d))
                
        Ire = {} # 
        for n in self.nodes:
            for d in demands.id:
                Ire[n, d] = model.addVar(vtype=GRB.BINARY, 
                   name='Ire_{}_{}'.format(n, d))
                
        III = {} # 
        for n in self.nodes:
            for d in demands.id:
                III[n, d] = model.addVar(vtype=GRB.BINARY, 
                   name='III_{}_{}'.format(n, d))
                
        I = {}
        for n in self.nodes:
            I[n] = model.addVar(vtype=GRB.BINARY, name='I_{}'.format(n))
                
        NNN = {}
        for n in self.nodes:
            NNN[n] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=Nmax, 
               name='NNN_{}'.format(n))
            
        X = {}
        for l in self.links.id:
            for d in demands.id:
                X[l, d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                 name='X_{}_{}'.format(l, d))
                
        Y = {}
        for n in self.nodes:
            for d in demands.id:
                Y[n, d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, 
                     ub=bigM1, name='Ynode_{}_{}'.format(n, d))
                
        UsageL = {} # if demand d uses link l
        for l in self.links.id:
            for d in demands.id:
                UsageL[l, d] = model.addVar(vtype=GRB.BINARY, 
                    name='UsageL_{}_{}'.format(l, d))
                
        Fstart = {} # the start frequency of demand d
        for d in demands.id:
            Fstart[d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                  name='Fstart_{}'.format(d))
            
        Delta = {} # order between demands
        for d1 in demands.id:
            for d2 in demands.id:
                Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                     name='Delta_{}_{}'.format(d1, d2))                    

        # GN variables  
        #    dvar float GNi[Links][Demands] in 0..10000;
        GNi={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                GNi[l,d]=model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=bigM1,
                         name='GNi_{}_{}'.format(l,d))
        
        #dvar float G1[Links][Demands] in 0..10000;
        G1={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                G1[l,d]=model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=bigM1,
                         name='G1_{}_{}'.format(l,d))

        #dvar float GASEws[Links] in 0..10000;
        GASEws={}
        for l in self.links.id:
            GASEws[l]=model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=bigM1,
                  name='GASEws_{}'.format(l))
             
        #dvar float GNliws[Links][Demands]in 0..10000;
        GNliws={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                GNliws[l,d]=model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=bigM1,
                      name='GNliws_{}_{}'.format(l,d))
                 
        #dvar float A1[Links][Demands] in 0..10000;
        A1={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                A1[l,d]=model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=bigM1,
                     name='A1_{}_{}'.format(l,d))
                 
        #var int UsageL1[Links][Demands][Demands]in 0..1;
        UsageL1={}# intermedia variables for product
        for l in self.links.id:
            for d1 in demands.id:
                for d2 in demands.id:
                    UsageL1[l,d1,d2]=model.addVar(vtype=GRB.BINARY,
                            name='UsaegL1_{}_{}_{}'.format(l,d1,d2))
                     
        #dvar float Asenli[Links][Demands]in 0..10000;
        Asenli={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                Asenli[l,d]=model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=bigM1,
                         name='Asenli_{}_{}'.format(l,d))
                 
        UseAsenli = {}
        for l in self.links.id:
            for d in demands.id:
                UseAsenli[l, d] = model.addVar(vtype=GRB.CONTINUOUS,lb=0,
                         ub=bigM1, name='UseAsenli_{}_{}'.format(l, d))
            
        Total = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=self.n_nodes, 
                             name='Total')
        
        c = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM2, name='c')
        
        model.update()
        
        # define constraints
        model.addConstr(Total==quicksum(I[n] for n in self.nodes), 
                        name='total')

        for d in demands.id:
            model.addConstr(c>=Fstart[d]+demands.data_rates[d], 
                            name='c_{}'.format(d))

        # flow conservation
        for n in self.nodes:
            for d in demands.id:
                model.addConstr(-quicksum(UsageL[l, d] for l in self.links.id
                                         if self.links.source[l]==n)+
                                quicksum(UsageL[l, d] for l in self.links.id
                                         if self.links.destination[l]==n)
                                == supply[n, d], 
                                name='flow_{}_{}'.format(n, d))

        for d1 in demands.id:
            for d2 in demands.id:
                if d1!=d2:
                    model.addConstr(Delta[d1, d2]+Delta[d2, d1]==1,
                                    name='delta_{}_{}'.format(d1, d2))
                else:
                    model.addConstr(Delta[d1, d2]+Delta[d2, d1]==0,
                                    name='delta_{}_{}'.format(d1, d2))
                    
        for d1 in demands.id:
            for d2 in demands.id:
                for l in self.links.id:
                    if d1!=d2:
                        model.addConstr(Fstart[d1]-Fstart[d2]<=
                                        bigM3*(3-Delta[d1, d2]-
                                               UsageL[l, d1]-UsageL[l, d2]),
                                        name='sharelink1_{}_{}_{}'.format(d1,d2,l))
                    
        for d1 in demands.id:
            for d2 in demands.id:
                for l in self.links.id:
                    if d1!=d2:
                        model.addConstr(Fstart[d1]-Fstart[d2]+
                                        demands.data_rates[d1]+G<=
                                        bigM3*(3-Delta[d1, d2]-UsageL[l, d1]-
                                        UsageL[l, d2]),
                                        name='sharelink2_{}_{}_{}'.format(d1,d2,l))
        
        for l in self.links.id:
            model.addConstr(GASEws[l]==self.links.length[l]*cofase,
                            name='gasews_{}'.format(l))
            
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(A1[l, d]<=bigM2*UsageL[l, d],
                                name='a11_{}_{}'.format(l,d))
                model.addConstr(A1[l, d]<=GASEws[l],name='a12_{}_{}'.format(l,d))
                model.addConstr(A1[l, d]>=GASEws[l]-(1-UsageL[l, d])*bigM2,
                                name='a13_{}_{}'.format(l,d))
                model.addConstr(A1[l, d]>=0,
                                name='a14_{}_{}'.format(l,d))
        
        for l in self.links.id:
            for d1 in demands.id:
                for d2 in demands.id:
                    model.addConstr(UsageL1[l, d1, d2]<=UsageL[l, d1],
                                    name='usagel11_{}_{}'.format(l,d1,d2))
                    model.addConstr(UsageL1[l, d1, d2]<=UsageL[l, d2],
                                    name='usagel12_{}_{}'.format(l,d1,d2))
                    model.addConstr(UsageL1[l, d1, d2]>=UsageL[l, d1]+
                                    UsageL[l, d2]-1,
                                    name='usagel13_{}_{}'.format(l,d1,d2))
                    
        for l in self.links.id:
            for d1 in demands.id:
                model.addConstr(GNi[l, d1]==np.log(rou*
                    demands.data_rates[d1]*demands.data_rates[d1])
#                       -UsageL[l, d1]*np.log((1.5*demands.data_rates[d1]+G)/
#                        (G+0.5*demands.data_rates[d1]))
                    +quicksum(UsageL1[l,d1,d3]*np.log(
                    (demands.data_rates[d3]+0.5*demands.data_rates[d1]+G)/
                    (G+0.5*demands.data_rates[d1]) ) 
                    for d3 in demands.id if d3!=d1),
                    name='usagel14_{}_{}'.format(l,d1))
        
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(GNliws[l, d]==self.links.length[l]*GNi[l, d],
                        name='gnliws_{}_{}'.format(l,d))
                
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(G1[l, d]<=bigM2*UsageL[l, d],
                        name='g11_{}_{}'.format(l,d))
                model.addConstr(G1[l, d]<=GNliws[l, d],
                        name='g12_{}_{}'.format(l,d))
                model.addConstr(G1[l, d]>=GNliws[l, d]-(1-UsageL[l, d])*bigM2,
                        name='g13_{}_{}'.format(l,d))
                model.addConstr(G1[l, d]>=0,
                        name='g14_{}_{}'.format(l,d))
                
        for n in self.nodes:
            for d in demands.id:
                model.addConstr(Y[n, d]<=self.Noise,
                        name='y_{}_{}'.format(n,d))
                
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(U[l, d]<=bigM2*UsageL[l, d],
                        name='u11_{}_{}'.format(l,d))
                model.addConstr(U[l, d]<=Y[self.links.source[l], d],
                        name='u12_{}_{}'.format(l,d))
                model.addConstr(U[l, d]>=Y[self.links.source[l], d]
                    -(1-UsageL[l, d])*bigM2,name='u13_{}_{}'.format(l,d))
                model.addConstr(U[l, d]>=0,name='u14_{}_{}'.format(l,d))
                
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(Asenli[l, d]==miu*G1[l, d]+A1[l, d],
                        name='asenli_{}_{}'.format(l,d))
                model.addConstr(UseAsenli[l, d]<=bigM1*UsageL[l, d],
                        name='useasenli11_{}_{}'.format(l,d))
                model.addConstr(UseAsenli[l, d]<=Asenli[l, d],
                        name='useasenli12_{}_{}'.format(l,d))
                model.addConstr(UseAsenli[l, d]>=Asenli[l, d]
                        -bigM1*(1-UsageL[l, d]),
                        name='useasenli13_{}_{}'.format(l,d))
                model.addConstr(UseAsenli[l, d]>=0,
                        name='useasenli14_{}_{}'.format(l,d))
                
        for n in self.nodes:
            for d in demands.id:
                model.addConstr(Y[n, d]==quicksum(X[l, d]+UseAsenli[l, d]
                    for l in self.links.id if self.links.destination[l]==n),
                    name='yx_{}_{}'.format(n, d))
                model.addConstr(III[n, d]==1-Ire[n, d], 
                                name='III_{}_{}'.format(n, d))
                
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(X[l, d]<=bigM1*III[self.links.source[l], d],
                        name='x11_{}_{}'.format(l,d))
                model.addConstr(X[l, d]<=U[l, d],
                        name='x12_{}_{}'.format(l,d))
                model.addConstr(X[l, d]>=U[l, d]-
                                bigM1*(1-III[self.links.source[l], d]),
                        name='x13_{}_{}'.format(l,d))
                model.addConstr(X[l, d]>=0,name='x14_{}_{}'.format(l,d))
        
        for n in self.nodes:
            model.addConstr(NNN[n]==quicksum(Ire[n, d] for d in demands.id),
                            name='nnn_{}'.format(n))
            model.addConstr(I[n]*Nmax>=NNN[n], name='nmax_{}'.format(n))
            
        # objective
        model.setObjective(c+epsilon_total*Total+epsilon_nnn*quicksum(NNN[n] 
            for n in self.nodes), GRB.MINIMIZE)

        # set gurobi parameters
#        if len(kwargs):
#            for key, value in kwargs.items():
#                try:
#                    setattr(model.params, key, value)
#                except:
#                    pass
        model.update()
        
        return model
        
    def create_model_tr(self, demands, **kwargs):
        '''Create MILP for TR at the first iteration of the stage, no demands 
        are fixed, do not optimize
        '''
        n_demands = demands.shape[0]
        
        # define supply 
        supply = np.zeros((self.n_nodes, n_demands))
        for n in range(self.n_nodes):
            for d in demands.id:
                if demands.source[d]==n:
                    supply[n, d] = -1

        for n in range(self.n_nodes):
            for d in demands.id:
                if demands.destination[d]==n:
                    supply[n, d] = 1
        
        model = Model('TR_{}'.format(demands.shape[0]))
        model.Params.UpdateMode = 1
        
        # define variables
        UsageL = {} # if demand d uses link l
        for l in self.links.id:
            for d in demands.id:
                UsageL[l, d] = model.addVar(vtype=GRB.BINARY, 
                    name='UsageL_{}_{}'.format(l, d))
                
        Fstart = {} # the start frequency of demand d
        for d in demands.id:
            Fstart[d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                  name='Fstart_{}'.format(d))
            
        Delta = {} # order between demands
        for d1 in demands.id:
            for d2 in demands.id:
                Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                     name='Delta_{}_{}'.format(d1, d2))
                    
        U = {} # U[a, b] = UsageL[a,b]*Ynode[a,b]
        for l in self.links.id:
            for d in demands.id:
                U[l, d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                 name='U_{}_{}'.format(l, d))
                
        Ire = {} # 
        for n in self.nodes:
            for d in demands.id:
                Ire[n, d] = model.addVar(vtype=GRB.BINARY, 
                   name='Ire_{}_{}'.format(n, d))
                
        III = {} # 
        for n in self.nodes:
            for d in demands.id:
                III[n, d] = model.addVar(vtype=GRB.BINARY, 
                   name='III_{}_{}'.format(n, d))
                
        I = {}
        for n in self.nodes:
            I[n] = model.addVar(vtype=GRB.BINARY, name='I_{}'.format(n))
                
        NNN = {}
        for n in self.nodes:
            NNN[n] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=Nmax, 
               name='NNN_{}'.format(n))
            
        X = {}
        for l in self.links.id:
            for d in demands.id:
                X[l, d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                 name='X_{}_{}'.format(l, d))
                
        Ynode = {}
        for n in self.nodes:
            for d in demands.id:
                Ynode[n, d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, 
                     ub=bigM1, name='Ynode_{}_{}'.format(n, d))
                
        Total = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=self.n_nodes, 
                             name='Total')
        
        c = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM2, name='c')
        
        model.update()
        
        # define constraints
        model.addConstr(Total==quicksum(I[n] for n in self.nodes))
        
        for d in demands.id:
            model.addConstr(c>=Fstart[d]+demands.data_rates[d])
            
        # flow conservation
        for n in self.nodes:
            for d in demands.id:
                model.addConstr(-quicksum(UsageL[l, d] for l in self.links.id
                                         if self.links.source[l]==n)+
                                quicksum(UsageL[l, d] for l in self.links.id
                                         if self.links.destination[l]==n)
                                == supply[n, d])
                                
        for d1 in demands.id:
            for d2 in demands.id:
                if d1!=d2:
                    model.addConstr(Delta[d1, d2]+Delta[d2, d1]==1)
                else:
                    model.addConstr(Delta[d1, d2]+Delta[d2, d1]==0)
                    
        for d1 in demands.id:
            for d2 in demands.id:
                for l in self.links.id:
                    if d1!=d2:
                        model.addConstr(Fstart[d1]-Fstart[d2]<=
                                        bigM3*(3-Delta[d1, d2]-
                                               UsageL[l, d1]-UsageL[l, d2]))
                    
        for d1 in demands.id:
            for d2 in demands.id:
                for l in self.links.id:
                    if d1!=d2:
                        model.addConstr(Fstart[d1]-Fstart[d2]+
                                        demands.data_rates[d1]+G<=
                                        bigM3*(3-Delta[d1, d2]-UsageL[l, d1]-
                                        UsageL[l, d2]))
                    
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(U[l, d]<=UsageL[l, d]*demands.TR[d])
                model.addConstr(U[l, d]<=Ynode[self.links.source[l], d])
                model.addConstr(Ynode[self.links.source[l], d]-U[l, d]<=
                                demands.TR[d]*(1-UsageL[l, d]))
                
        for n in self.nodes:
            for d in demands.id:
                model.addConstr(Ynode[n, d]==
                                quicksum(X[l, d]+UsageL[l, d]*
                                          self.links.length[l] 
                                          for l in self.links.id 
                                          if self.links.destination[l]==n))
                
        for n in self.nodes:
            for d in demands.id:
                model.addConstr(III[n, d]==1-Ire[n, d])
                
        for l in self.links.id:
            for d in demands.id:
                model.addConstr(X[l, d]<=bigM1*III[self.links.source[l], d])
                model.addConstr(X[l, d]<=U[l, d])
                model.addConstr(X[l, d]>=U[l, d]-
                                bigM1*(1-III[self.links.source[l], d]))
                model.addConstr(X[l, d]>=0)
        
        for n in self.nodes:
            model.addConstr(NNN[n]==quicksum(Ire[n, d] for d in demands.id))
            model.addConstr(I[n]*Nmax>=NNN[n])
            
        # objective
        model.setObjective(c+epsilon_total*Total+epsilon_nnn*quicksum(NNN[n] 
            for n in self.nodes), GRB.MINIMIZE)
            
        # set gurobi parameters
#        if len(kwargs):
#            for key, value in kwargs.items():
#                try:
#                    setattr(model.params, key, value)
#                except:
#                    pass
        model.update()

        return model
    
    def solve_model_all(self, model, demands, num_resolve=1, **kwargs):
        '''Solve MILP model with all demands as variables, used for the first
        iteration at the first stage
        '''
        if len(kwargs):
            for key, value in kwargs.items():
                try:
                    setattr(model.params, key, value)
                except:
                    pass
        
        model.optimize()

        # save files
        if len(kwargs):
            for key, value in kwargs.items():
                if key=='write':
                    if type(value) is list:
                        for i in value:
                            try:
                                model.write(i)
                            except:
                                pass
                    elif type(value) is str:
                        try:
                            model.write(value)
                        except:
                            pass
              
        if model.SolCount>=1:
            flag_success = True
            # construct solutions for UsageL and Delta
            UsageLx = {}
            for l in self.links.id:
                for d in demands.id:
                    UsageL_ = model.getVarByName('UsageL_{}_{}'.format(l, d))
                    try:
                        if UsageL_.x<0.5:
                            UsageLx[l, d] = 0
                        else:
                            UsageLx[l, d] = 1
                    except:
                        print(l, d)
                        raise
                           
            Deltax = {}
            for d1 in demands.id:
                for d2 in demands.id:
                    Delta_ = model.getVarByName('Delta_{}_{}'.format(d1, d2))
                    if Delta_.x <0.5:
                        Deltax[d1, d2] = 0
                    else:
                        Deltax[d1, d2] = 1
                              
            Fstartx = {}
            for d in demands.id:
                Fstart_ = model.getVarByName('Fstart_{}'.format(d))
                Fstartx[d] = np.round(Fstart_.x*2)/2
                          
            Totalx = round(model.getVarByName('Total').x, 2)
            cx = round(model.getVarByName('c').x, 2)
                                   
            NNNx = {}
            for n in self.nodes:
                NNN_ = model.getVarByName('NNN_{}'.format(n))
                NNNx[n] = round(NNN_.x, 2)
        else:
            flag_success = False
            # solving is failed
            UsageLx = {(l, d):0 for l in self.links.id for d in demands.id}
            Deltax = {(d1, d2):0 for d1 in demands.id for d2 in demands.id}                              
            Fstartx = {d:0 for d in demands.id}  
            Totalx = INF
            cx = INF
            NNNx = {n:0 for n in self.nodes}
        
        solutions = {}
        solutions['UsageL'] = UsageLx
        solutions['Fstart'] = Fstartx
        solutions['Delta'] = Deltax
        solutions['Total'] = Totalx
        solutions['c'] = cx
        solutions['flag_success'] = flag_success
        solutions['NNN'] = NNNx
        
        return model, solutions
        
    
    def modify_model(self, model, demands, previous_solutions, 
                     num_resolve=1, miphint=True, **kwargs):
        '''Modify MILP model and solve, VarName0 gives hint while VarName gives
        lb and ub
        '''
        demands_added = previous_solutions['demands_added'] # new demands
        demands_fixed = previous_solutions['demands_fixed'] # old demands
        UsageLx0 = previous_solutions['UsageL0']
        Deltax0 = previous_solutions['Delta0']
        Fstartx0 = previous_solutions['Fstart0']
        UsageLx = previous_solutions['UsageL']
        Deltax = previous_solutions['Delta']
        ObjVal = round(100*previous_solutions['ObjVal'], 2)/100
        flag_success = previous_solutions['flag_success']
        demands_all = list(set(demands_added).union(set(demands_fixed)))
        demands = demands.loc[demands.id.isin(demands_all), :]

        # modify variable bounds
        for l in self.links.id:
            for d in demands.id:
                UsageL_ = model.getVarByName('UsageL_{}_{}'.format(l, d))
                if ((d in demands_fixed) and ((l, d) in UsageLx.keys()) 
                    and flag_success):
                    UsageL_.lb = UsageLx[l, d]
                    UsageL_.ub = UsageLx[l, d]
                else:
                    UsageL_.lb = 0
                    UsageL_.ub = 1
                    if miphint and ((l, d) in UsageLx0.keys()):
                        UsageL_.VarHintVal = UsageLx0[l, d]
                    
        for d in demands.id:
            Fstart_ = model.getVarByName('Fstart_{}'.format(d))
            if miphint and (d in Fstartx0.keys()):
                Fstart_.VarHintVal = Fstartx0[d]
                
        for d1 in demands.id:
            for d2 in demands.id:
                Delta_ = model.getVarByName('Delta_{}_{}'.format(d1, d2))
                if ((d1!=d2) and (d1 in demands_fixed) 
                    and (d2 in demands_fixed) and flag_success):
                    Delta_.lb = Deltax[d1, d2]
                    Delta_.ub = Deltax[d1, d2]
                elif d1==d2:
                    Delta_.lb = 0
                    Delta_.ub = 0
                elif d1!=d2:
                    Delta_.lb = 0
                    Delta_.ub = 1                
                    if miphint and ((d1, d2) in Deltax0.keys()):
                        Delta_.VarHintVal = Deltax0[d1, d2]
                
        # set gurobi parameters
        if len(kwargs):
            for key, value in kwargs.items():
                try:
                    setattr(model.params, key, value)
                except:
                    pass
        model.update()

        model.optimize()
        while ((model.SolCount<1 or round(100*model.ObjVal, 2)/100>ObjVal) 
            and num_resolve)>0:
            if len(kwargs):
                for key, value in kwargs.items():
                    try:
                        setattr(model.params, key, value)
                    except:
                        pass
            model.update()
            model.optimize()
            num_resolve -= 1

        # save files
        if len(kwargs):
            for key, value in kwargs.items():
                if key=='write':
                    if type(value) is list:
                        for i in value:
                            try:
                                model.write(i)
                            except:
                                pass
                    elif type(value) is str:
                        try:
                            model.write(value)
                        except:
                            pass
              
        if model.SolCount>=1:
            flag_success = True
            # construct solutions for UsageL and Delta
            UsageLx = {}
            for l in self.links.id:
                for d in demands.id:
                    UsageL_ = model.getVarByName('UsageL_{}_{}'.format(l, d))
                    if UsageL_.x<0.5:
                        UsageLx[l, d] = 0
                    else:
                        UsageLx[l, d] = 1
                           
            Deltax = {}
            for d1 in demands.id:
                for d2 in demands.id:
                    Delta_ = model.getVarByName('Delta_{}_{}'.format(d1, d2))
                    if Delta_.x <0.5:
                        Deltax[d1, d2] = 0
                    else:
                        Deltax[d1, d2] = 1
                              
            Fstartx = {}
            for d in demands.id:
                Fstart_ = model.getVarByName('Fstart_{}'.format(d))
                Fstartx[d] = np.round(Fstart_.x*2)/2
                          
            Totalx = round(model.getVarByName('Total').x, 2)
            cx = round(model.getVarByName('c').x, 2)
                                   
            NNNx = {}
            for n in self.nodes:
                NNN_ = model.getVarByName('NNN_{}'.format(n))
                NNNx[n] = round(NNN_.x, 2)
                    
        else:
            flag_success = False
            # solving is failed
            UsageLx = {(l, d):0 for l in self.links.id for d in demands.id}
            Deltax = {(d1, d2):0 for d1 in demands.id for d2 in demands.id}                              
            Fstartx = {d:0 for d in demands.id}  
            Totalx = INF
            cx = INF
            NNNx = {n:0 for n in self.nodes}
        
        solutions = {}
        solutions['UsageL'] = UsageLx
        solutions['Fstart'] = Fstartx
        solutions['Delta'] = Deltax
        solutions['Total'] = Totalx
        solutions['c'] = cx
        solutions['NNN'] = NNNx
        solutions['flag_success'] = flag_success
        
        return model, solutions
    
    def scheduler(self, demands, idx_iter, idx_stage, n_demands_initial, 
                  max_added, previous=None):
        '''Randomly hold out some demands for re-optimization
        '''
        if (idx_iter==0) and (idx_stage==0):
            # the first iteration in the first stage
            demands_added = demands.iloc[:n_demands_initial].id.values.tolist()
            demands_fixed = []
        elif (idx_iter==0) and (idx_stage>=1):
            # the first iteration of a non-first stage
            n_demands_pre = n_demands_initial+(idx_stage-1)*n_demands_per_stage
            n_demands_all = n_demands_initial+idx_stage*n_demands_per_stage
            demands_added = demands.iloc[n_demands_pre:
                n_demands_all].id.values.tolist()
            demands_fixed = demands.iloc[:n_demands_pre].id.values.tolist()
        elif idx_iter>=1:
            # the non-first iteration of a stage
            n_demands_all = n_demands_initial+idx_stage*n_demands_per_stage
            demands_all_id = demands.iloc[:n_demands_all].id.values.tolist()
            np.random.shuffle(demands_all_id)
            if previous is None:
                n_added = np.random.randint(1, max_added+1)
            elif previous['better']:
                n_added = np.random.randint(previous['n_added'], 
                    previous['n_added']+3)
            else:
                n_added = np.random.randint(previous['n_added']-2, 
                    previous['n_added']+1)
            n_demands_pre = n_demands_all-n_added
            demands_fixed = demands_all_id[:n_demands_pre]
            demands_added = demands_all_id[n_demands_pre:n_demands_all]
        
        idx = idx_stage*n_iter_per_stage+idx_iter
        timelimit = max(timelimit_baseline, 
            timelimit0*2**(np.floor(idx/n_iter_per_stage)/time_factor))
        
        return demands_added, demands_fixed, timelimit

    def iterate(self, demands, max_added=n_demands_per_stage, miphint=True, 
                num_resolve=[1, 1], **kwargs):
        '''Solve TR and GN simultaneously, the best solution from TR and GN
            is input to the next iteration of solving (either TR or GN) 
            as the start point
        '''
        tic = time.clock()
        
        n_demands = demands.shape[0]
        n_stages = int(np.ceil(n_demands/n_demands_per_stage))
        n_demands_initial = n_demands-(n_stages-1)*n_demands_per_stage
        num_resolve_gn = num_resolve[0]
        num_resolve_tr = num_resolve[1]
        
        iteration_history_tr = {}
        iteration_history_gn = {}
        idx = 0
        logfile = 'gurobi.log'
        for k in kwargs.keys():
            if k.lower()=='logfile':
                logfile = kwargs[k]
                break
        
        for idx_stage in range(n_stages):
            n_demands_in_stage = n_demands_initial+\
                idx_stage*n_demands_per_stage
            demands_in_stage = demands.iloc[:n_demands_in_stage, :]
            model_gn = self.create_model_gn(demands_in_stage, **kwargs)
            model_tr = self.create_model_tr(demands_in_stage, **kwargs)
            
            for idx_iter in range(n_iter_per_stage):
                if callable(max_added):
                    max_added_ = max_added(n_demands_in_stage)
                else:
                    max_added_ = max_added

                if idx<2:
                    previous = None
                else:
                    previous = {}
                    previous['n_added'] = len(iteration_history_gn[idx-1][
                        'demands_added'])
                    if (iteration_history_gn[idx-1]['ObjVal']<=
                        iteration_history_gn[idx-2]['ObjVal']):
                        previous['better'] = True
                    else:
                        previous['better'] = False

                demands_added, demands_fixed, timelimit = \
                    self.scheduler(demands, idx_iter, idx_stage, 
                    n_demands_initial, max_added_, previous=previous)

                if (idx_iter==0) and (idx_stage==0):
                    # the first iteration in the first stage
                    # solve TR
                    with open(logfile, 'a') as f:
                        f.write('\n#######################################\n')
                        f.write('TR: iteration {} at stage {}\n'.format(idx_iter, idx_stage))
                        f.write('Added demands: {}\n'.format(demands_added))
                        f.write('Fixed demands: {}\n'.format(demands_fixed))
                    
                    model_tr, solutions_tr = self.solve_model_all(model_tr,
                        demands_in_stage, num_resolve=num_resolve_tr, 
                        timelimit=timelimit, **kwargs)
                    toc_now = time.clock()
                    iteration_history_tr[idx] = {}
                    iteration_history_tr[idx]['elapsed_time'] = toc_now-tic
                    iteration_history_tr[idx]['demands_fixed'] = []
                    iteration_history_tr[idx]['demands_added'] = \
                        demands_in_stage.id.values.tolist()
                    iteration_history_tr[idx]['demands_solved'] = \
                        demands_in_stage.id.values.tolist()
                    iteration_history_tr[idx]['solutions'] = solutions_tr
                    iteration_history_tr[idx]['ObjVal'] = \
                        round(model_tr.ObjVal*100, 2)/100
                    iteration_history_tr[idx]['flag_success'] = \
                        solutions_tr['flag_success']
                        
                    # solve GN
                    with open(logfile, 'a') as f:
                        f.write('\n#######################################\n')
                        f.write('GN: iteration {} at stage {}\n'.format(
                            idx_iter, idx_stage))
                        f.write('Added demands: {}\n'.format(demands_added))
                        f.write('Fixed demands: {}\n'.format(demands_fixed))
                    
                    model_gn, solutions_gn = self.solve_model_all(model_gn, 
                        demands_in_stage, num_resolve=num_resolve_gn, 
                        timelimit=timelimit, **kwargs)
                    toc_now = time.clock()
                    iteration_history_gn[idx] = {}
                    iteration_history_gn[idx]['elapsed_time'] = toc_now-tic
                    iteration_history_gn[idx]['demands_fixed'] = []
                    iteration_history_gn[idx]['demands_added'] = \
                        demands_in_stage.id.values.tolist()
                    iteration_history_gn[idx]['demands_solved'] = \
                        demands_in_stage.id.values.tolist()
                    iteration_history_gn[idx]['solutions'] = solutions_gn
                    iteration_history_gn[idx]['ObjVal'] = \
                        round(model_gn.ObjVal*100, 2)/100
                    iteration_history_gn[idx]['flag_success'] = \
                        solutions_gn['flag_success']
                                        
                    idx += 1
                                        
                    continue
                
                # solve TR
                with open(logfile, 'a') as f:
                    f.write('\n#######################################\n')
                    f.write('TR: iteration {} at stage {}\n'.format(
                        idx_iter, idx_stage))
                    f.write('Added demands: {}\n'.format(demands_added))
                    f.write('Fixed demands: {}\n'.format(demands_fixed))
                
                previous_solutions_tr = {}
                previous_solutions_tr['demands_added'] = demands_added
                previous_solutions_tr['demands_fixed'] = demands_fixed
                previous_solutions_tr['UsageL'] = \
                    iteration_history_tr[idx-1]['solutions']['UsageL']
                previous_solutions_tr['Delta'] = \
                    iteration_history_tr[idx-1]['solutions']['Delta']
                previous_solutions_tr['ObjVal'] = \
                    iteration_history_tr[idx-1]['ObjVal']
                previous_solutions_tr['flag_success'] = \
                    iteration_history_tr[idx-1]['solutions']['flag_success']
                if (iteration_history_tr[idx-1]['ObjVal'] 
                    <= iteration_history_gn[idx-1]['ObjVal']):
                    previous_solutions_tr['UsageL0'] = \
                        iteration_history_tr[idx-1]['solutions']['UsageL']
                    previous_solutions_tr['Delta0'] = \
                        iteration_history_tr[idx-1]['solutions']['Delta']
                    previous_solutions_tr['Fstart0'] = \
                        iteration_history_tr[idx-1]['solutions']['Fstart']
                else:
                    previous_solutions_tr['UsageL0'] = \
                        iteration_history_gn[idx-1]['solutions']['UsageL']
                    previous_solutions_tr['Delta0'] = \
                        iteration_history_gn[idx-1]['solutions']['Delta']
                    previous_solutions_tr['Fstart0'] = \
                        iteration_history_gn[idx-1]['solutions']['Fstart']
                        
                if not previous_solutions_tr['flag_success']:
                    timelimit *= idx_stage+1
                model_tr, solutions_tr = self.modify_model(model_tr, 
                    demands, previous_solutions_tr, num_resolve=num_resolve_tr,
                    miphint=miphint, timelimit=timelimit, **kwargs)
                toc_now = time.clock()
                iteration_history_tr[idx] = {}
                iteration_history_tr[idx]['elapsed_time'] = toc_now-tic
                iteration_history_tr[idx]['demands_fixed'] = demands_fixed
                iteration_history_tr[idx]['demands_added'] = demands_added
                iteration_history_tr[idx]['demands_solved'] = \
                    list(set(demands_fixed).union(set(demands_added)))
                iteration_history_tr[idx]['solutions'] = solutions_tr
                iteration_history_tr[idx]['ObjVal'] = \
                    round(model_tr.ObjVal*100, 2)/100
                    
                # solve GN
                with open(logfile, 'a') as f:
                    f.write('\n#######################################\n')
                    f.write('GN: iteration {} at stage {}\n'.format(
                            idx_iter, idx_stage))
                    f.write('Added demands: {}\n'.format(demands_added))
                    f.write('Fixed demands: {}\n'.format(demands_fixed))
                
                previous_solutions_gn = {}
                previous_solutions_gn['demands_added'] = demands_added
                previous_solutions_gn['demands_fixed'] = demands_fixed
                previous_solutions_gn['UsageL'] = \
                    iteration_history_gn[idx-1]['solutions']['UsageL']
                previous_solutions_gn['Delta'] = \
                    iteration_history_gn[idx-1]['solutions']['Delta']
                previous_solutions_gn['ObjVal'] = \
                    iteration_history_gn[idx-1]['ObjVal']
                previous_solutions_gn['flag_success'] = \
                    iteration_history_gn[idx-1]['solutions']['flag_success']
                if (iteration_history_gn[idx-1]['ObjVal'] 
                    <= iteration_history_tr[idx]['ObjVal']):
                    previous_solutions_gn['UsageL0'] = \
                        iteration_history_gn[idx-1]['solutions']['UsageL']
                    previous_solutions_gn['Delta0'] = \
                        iteration_history_gn[idx-1]['solutions']['Delta']
                    previous_solutions_gn['Fstart0'] = \
                        iteration_history_gn[idx-1]['solutions']['Fstart']
                else:
                    previous_solutions_gn['UsageL0'] = \
                        iteration_history_tr[idx]['solutions']['UsageL']
                    previous_solutions_gn['Delta0'] = \
                        iteration_history_tr[idx]['solutions']['Delta']
                    previous_solutions_gn['Fstart0'] = \
                        iteration_history_tr[idx]['solutions']['Fstart']
                        
                if not previous_solutions_gn['flag_success']:
                    timelimit *= idx_stage+1
                model_gn, solutions_gn = self.modify_model(model_gn, 
                    demands, previous_solutions_gn, num_resolve=num_resolve_gn,
                    miphint=miphint, timelimit=timelimit, **kwargs)
                toc_now = time.clock()
                iteration_history_gn[idx] = {}
                iteration_history_gn[idx]['elapsed_time'] = toc_now-tic
                iteration_history_gn[idx]['demands_fixed'] = demands_fixed
                iteration_history_gn[idx]['demands_added'] = demands_added
                iteration_history_gn[idx]['demands_solved'] = \
                    list(set(demands_fixed).union(set(demands_added)))
                iteration_history_gn[idx]['solutions'] = solutions_gn
                iteration_history_gn[idx]['ObjVal'] = \
                    round(model_gn.ObjVal*100, 2)/100
                
                idx += 1
                
                # cheating...
                if (iteration_history_gn[idx-1]['ObjVal']>
                    iteration_history_tr[idx-1]['ObjVal']+0.5):
                    with open(logfile, 'a') as f:
                        f.write('\n#######################################\n')
                        f.write('I AM CHEATING!!! at stage {}\n'.format(
                            idx_stage))
                    previous_solutions_gn['demands_added'] = []
                    previous_solutions_gn['demands_fixed'] = \
                        list(set(demands_fixed).union(set(demands_added)))
                    previous_solutions_gn['UsageL'] = \
                        iteration_history_tr[idx-1]['solutions']['UsageL']
                    previous_solutions_gn['Delta'] = \
                        iteration_history_tr[idx-1]['solutions']['Delta']
                    previous_solutions_gn['ObjVal'] = \
                        iteration_history_tr[idx-1]['ObjVal']
                    previous_solutions_gn['flag_success'] = \
                        iteration_history_tr[idx-1]['solutions'][
                        'flag_success']
                    previous_solutions_gn['UsageL0'] = \
                        iteration_history_tr[idx-1]['solutions']['UsageL']
                    previous_solutions_gn['Delta0'] = \
                        iteration_history_tr[idx-1]['solutions']['Delta']
                    previous_solutions_gn['Fstart0'] = \
                        iteration_history_tr[idx-1]['solutions']['Fstart']
                            
                    if not previous_solutions_gn['flag_success']:
                        timelimit *= idx_stage+1
                        
                    model_gn, solutions_gn = self.modify_model(model_gn, 
                        demands, previous_solutions_gn, 
                        num_resolve=num_resolve_gn,
                        miphint=miphint, timelimit=timelimit, **kwargs)
                    toc_now = time.clock()
                    iteration_history_gn[idx] = {}
                    iteration_history_gn[idx]['elapsed_time'] = toc_now-tic
                    iteration_history_gn[idx]['demands_fixed'] = demands_fixed
                    iteration_history_gn[idx]['demands_added'] = demands_added
                    iteration_history_gn[idx]['demands_solved'] = \
                        list(set(demands_fixed).union(set(demands_added)))
                    iteration_history_gn[idx]['solutions'] = solutions_gn
                    iteration_history_gn[idx]['ObjVal'] = model_gn.ObjVal
    
                    iteration_history_tr[idx] = {}
                    iteration_history_tr[idx]['elapsed_time'] = \
                        iteration_history_tr[idx-1]['elapsed_time']
                    iteration_history_tr[idx]['demands_fixed'] = \
                        iteration_history_tr[idx-1]['demands_fixed']
                    iteration_history_tr[idx]['demands_added'] = \
                        iteration_history_tr[idx-1]['demands_added']
                    iteration_history_tr[idx]['demands_solved'] = \
                        iteration_history_tr[idx-1]['demands_solved']
                    iteration_history_tr[idx]['solutions'] = \
                        iteration_history_tr[idx-1]['solutions']
                    iteration_history_tr[idx]['ObjVal'] = \
                        iteration_history_tr[idx-1]['ObjVal']
                        
                    idx += 1

        toc = time.clock()
        self.total_runtime = toc-tic

        return iteration_history_tr, iteration_history_gn

def extract_history(iteration_history, variable_name):
    '''Extract the history of a certain variable in the iteration_history
    '''
    if variable_name in iteration_history[0]['solutions']:
        var_history = [iteration_history[i]['solutions'].get(variable_name)
            for i in range(len(iteration_history))]
    elif variable_name in iteration_history[0]:
        var_history = [iteration_history[i].get(variable_name)
            for i in range(len(iteration_history))]
#    elif hasattr(iteration_history[0]['model'], variable_name):
#        var_history = [getattr(iteration_history[i]['model'], variable_name)
#            for i in range(len(iteration_history))]

    return var_history

def read_demands(demands_csv, modulation='bpsk'):
    '''read demand csv file'''
    demands = pd.read_csv(demands_csv, header=None)
    demands.reset_index(drop=False, inplace=True)
    demands.columns = ['id', 'source', 'destination', 'data_rates']
    n_demands = demands.shape[0]
    # choose modulation format
    if modulation=='qpsk':
        bpsk_tr = pd.read_csv('qpsk_TR.csv', header=None)
        bpsk_tr.columns = ['data_rate', 'distance']
        bpsk_tr.distance = bpsk_tr.distance/100
        bpsk_tr.set_index('data_rate', inplace=True)
        tr = [float(bpsk_tr.loc[int(np.round(demands.data_rates[i]))]) 
            for i in range(n_demands)]
    elif modulation=='bpsk':
        bpsk_tr = pd.read_csv('bpsk_TR.csv', header=None)
        bpsk_tr.columns = ['data_rate', 'distance']
        bpsk_tr.distance = bpsk_tr.distance/100
        bpsk_tr.set_index('data_rate', inplace=True)
        tr = [float(bpsk_tr.loc[int(np.round(demands.data_rates[i]))]) 
            for i in range(n_demands)]
        
    demands['TR'] = tr

    return demands

def save_data(file_name, data):
    """File name must ends with .josn
    """
    with open(file_name, 'wb') as f:
        pickle.dump(data, f)
        
def read_data(file_name):
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
        
    return data

def half_n(n):
    return int(round(n/2))

def whole_n(n):
    return n