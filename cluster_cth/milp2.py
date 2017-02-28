# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 17:35:05 2017

@author: yx4vf
"""

import numpy as np
import pandas as pd
from gurobipy import *
import time
import copy
import pickle
import gc

# fiber parameters
INF = np.inf # infinity
G = 12.5 # guardband
Nmax = 10 # max number of regenerator circuits per regenerator node
cofase = 23.86 # ASE coefficient
rou = 2.11*10**-3
miu = 1.705


# modelling parameters
bigM1 = 10**5 
bigM2 = 10**5
bigM3 = 2*10**6 

# scheduler parameters
n_demands_initial = 5
n_iter_per_stage = 6 # 10
th_mipgap = 0.01
n_demands_increment = 5
timelimit_baseline = 960
timelimit0 = 60
time_factor = 1.3

np.random.seed(0) # set random seed

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
    
    def solve_all_gn(self, demands, **kwargs):
        '''Formulate and solve
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
        
        tic = time.clock()
        model = Model('TR')
        
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
            NNN[n] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=10, 
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
        model.setObjective(c+Total, GRB.MINIMIZE)
            
        # set gurobi parameters
        if len(kwargs):
            for key, value in kwargs.items():
                try:
                    setattr(model.params, key, value)
                except:
                    pass
        model.update()
                
        toc = time.clock()
        self.model_time = toc-tic
        
        model.optimize()
        
        toc2 = time.clock()
        self.solve_time = toc2-toc
        
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
                        
        try:
            # construct solutions for UsageL and Delta
            UsageLx = {}
            for l in self.links.id:
                for d in demands.id:
                    if UsageL[l, d].x<0.5:
                        UsageLx[l, d] = 0
                    else:
                        UsageLx[l, d] = 1
                           
            Deltax = {}
            for d1 in demands.id:
                for d2 in demands.id:
                    if Delta[d1, d2].x <0.5:
                        Deltax[d1, d2] = 0
                    else:
                        Deltax[d1, d2] = 1
    
            Fstartx = {}
            for d in demands.id:
                Fstartx[d] = np.round(Fstart[d].x*2)/2
                       
            Ux = {} # U[a, b] = UsageL[a,b]*Ynode[a,b]
            for l in self.links.id:
                for d in demands.id:
                    Ux[l, d] = U[l, d].x
                      
            Irex = {} # 
            for n in self.nodes:
                for d in demands.id:
                    Irex[n, d] = Ire[n, d].x
                        
            IIIx = {} # 
            for n in self.nodes:
                for d in demands.id:
                    IIIx[n, d] = III[n, d].x
                        
            Ix = {}
            for n in self.nodes:
                Ix[n] = I[n].x
                  
            NNNx = {}
            for n in self.nodes:
                NNNx[n] = NNN[n].x
                    
            Xx = {}
            for l in self.links.id:
                for d in demands.id:
                    Xx[l, d] = X[l, d].x
                      
            Yx = {}
            for n in self.nodes:
                for d in demands.id:
                    Yx[n, d] = Y[n, d].x
                          
            Totalx = Total.x
            
            cx = c.x
            
            GNix = {}
            for l in self.links.id:
                for d in demands.id:
                    GNix[l, d] = GNi[l, d].x
                        
            G1x={}# intermedia variables for product
            for l in self.links.id:
                for d in demands.id:
                    G1x[l,d] = G1[l,d].x
    
            #dvar float GASEws[Links] in 0..10000;
            GASEwsx={}
            for l in self.links.id:
                GASEwsx[l] = GASEws[l].x
                 
            #dvar float GNliws[Links][Demands]in 0..10000;
            GNliwsx={}# intermedia variables for product
            for l in self.links.id:
                for d in demands.id:
                    GNliwsx[l,d] = GNliws[l,d].x
                     
            #dvar float A1[Links][Demands] in 0..10000;
            A1x={}# intermedia variables for product
            for l in self.links.id:
                for d in demands.id:
                    A1x[l,d] = A1[l,d].x
                     
            #var int UsageL1[Links][Demands][Demands]in 0..1;
            UsageL1x={}# intermedia variables for product
            for l in self.links.id:
                for d1 in demands.id:
                    for d2 in demands.id:
                        UsageL1x[l,d1,d2] = UsageL1[l,d1,d2].x
                         
            #dvar float Asenli[Links][Demands]in 0..10000;
            Asenlix={}# intermedia variables for product
            for l in self.links.id:
                for d in demands.id:
                    Asenlix[l,d] = Asenli[l,d].x
                     
            UseAsenlix = {}
            for l in self.links.id:
                for d in demands.id:
                    UseAsenlix[l, d] = UseAsenli[l, d].x
                              
            G1x={}# intermedia variables for product
            for l in self.links.id:
                for d in demands.id:
                    G1x[l,d]=G1[l, d].x
            
            solutions = {}
            solutions['UsageL'] = UsageLx
            solutions['Fstart'] = Fstartx
            solutions['Delta'] = Deltax
            solutions['U'] = Ux
            solutions['Ire'] = Irex
            solutions['III'] = IIIx
            solutions['I'] = Ix
            solutions['NNN'] = NNNx
            solutions['X'] = Xx
            solutions['Y'] = Yx
            solutions['Total'] = Totalx
            solutions['c'] = cx
            solutions['GASEws'] = GASEwsx
            solutions['GNliws'] = GNliwsx
            solutions['A1'] = A1x
            solutions['UsageL1'] = UsageL1x
            solutions['Asenli'] = Asenlix
            solutions['UseAsenli'] = UseAsenlix
            solutions['G1'] = G1x
            
            return model, solutions, UsageLx, Deltax
        
        except:
            return model
        
    def solve_all_tr(self, demands, **kwargs):
        '''Formulate and solve
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
        
        tic = time.clock()
        model = Model('TR')
        
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
            NNN[n] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=10, 
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
        model.setObjective(c+Total, GRB.MINIMIZE)
            
        # set gurobi parameters
        if len(kwargs):
            for key, value in kwargs.items():
                try:
                    setattr(model.params, key, value)
                except:
                    pass
        model.update()
                
        toc = time.clock()
        self.model_time = toc-tic
        
        model.optimize()
        
        toc2 = time.clock()
        self.solve_time = toc2-toc
        
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
                        
        try:
        # construct solutions for UsageL and Delta
            UsageLx = {}
            for l in self.links.id:
                for d in demands.id:
                    if UsageL[l, d].x<0.5:
                        UsageLx[l, d] = 0
                    else:
                        UsageLx[l, d] = 1
                           
            Deltax = {}
            for d1 in demands.id:
                for d2 in demands.id:
                    if Delta[d1, d2].x <0.5:
                        Deltax[d1, d2] = 0
                    else:
                        Deltax[d1, d2] = 1
    
            Fstartx = {}
            for d in demands.id:
                Fstartx[d] = np.round(Fstart[d].x*2)/2
                       
            Ux = {} # U[a, b] = UsageL[a,b]*Ynode[a,b]
            for l in self.links.id:
                for d in demands.id:
                    Ux[l, d] = U[l, d].x
                      
            Irex = {} # 
            for n in self.nodes:
                for d in demands.id:
                    Irex[n, d] = Ire[n, d].x
                        
            IIIx = {} # 
            for n in self.nodes:
                for d in demands.id:
                    IIIx[n, d] = III[n, d].x
                        
            Ix = {}
            for n in self.nodes:
                Ix[n] = I[n].x
                  
            NNNx = {}
            for n in self.nodes:
                NNNx[n] = NNN[n].x
                    
            Xx = {}
            for l in self.links.id:
                for d in demands.id:
                    Xx[l, d] = X[l, d].x
                      
            Ynodex = {}
            for n in self.nodes:
                for d in demands.id:
                    Ynodex[n, d] = Ynode[n, d].x
                          
            Totalx = Total.x
            
            cx = c.x
            
            solutions = {}
            solutions['UsageL'] = UsageLx
            solutions['Fstart'] = Fstartx
            solutions['Delta'] = Deltax
            solutions['U'] = Ux
            solutions['Ire'] = Irex
            solutions['III'] = IIIx
            solutions['I'] = Ix
            solutions['NNN'] = NNNx
            solutions['X'] = Xx
            solutions['Ynode'] = Ynodex
            solutions['Total'] = Totalx
            solutions['c'] = cx
        
            return model, solutions, UsageLx, Deltax
        
        except:
            return model
        
    def solve_partial_gn(self, demands, previous_solutions, mipstart=False, 
                      num_solve=2, FeasibilityTol=1e-9, IntFeasTol=1e-9, 
                      OptimalityTol=1e-9, **kwargs):
        '''Formulate and solve iteratively
        previous_solutions is dict, contains:
            - UsageL from the previous solve, dict
            - Delta from the previous solve, dict
            - demands_added new demands that will be allocated this time, list
            - demands_fixed old demands allocated previously, list, 
                can be empty
            we should make sure that UsageL and Delta contain solutions for 
            all the demands in demands_fixed
        '''
        # process input data
        demands_added = previous_solutions['demands_added'] # new demands
        demands_fixed = previous_solutions['demands_fixed'] # old demands
        UsageLx0 = previous_solutions['UsageL0']
        Deltax0 = previous_solutions['Delta0']
        Fstartx0 = previous_solutions['Fstart0']
        UsageLx = previous_solutions['UsageL']
        Deltax = previous_solutions['Delta']
        ObjVal = previous_solutions['ObjVal']
        demands_all = list(set(demands_added).union(set(demands_fixed)))
#        print(demands_added)
#        print(demands_fixed)
#        print(demands_all)
        n_demands = len(demands_all)
        demands = demands.loc[demands.id.isin(demands_all), :]
        

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
                          
        
        tic = time.clock()
        model = Model('TR')
        
        # define variables
        UsageL = {} # if demand d uses link l
        for l in self.links.id:
            for d in demands.id:
                if d in demands_fixed:
                    UsageL[l, d] = model.addVar(vtype=GRB.BINARY, 
                        name='UsageL_{}_{}'.format(l, d), 
                        lb=UsageLx[l, d], ub=UsageLx[l, d])
                elif d in demands_added:
                    UsageL[l, d] = model.addVar(vtype=GRB.BINARY, 
                        name='UsageL_{}_{}'.format(l, d))
                if mipstart:
                    try:
                        UsageL[l, d].start=UsageLx0[l, d]
                    except:
                        pass
                
        Fstart = {} # the start frequency of demand d
        for d in demands.id:
            Fstart[d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                  name='Fstart_{}'.format(d))
            if mipstart:
                try:
                    Fstart[d].start = Fstartx0[d]
                except:
                    pass
            
        Delta = {} # order between demands
        for d1 in demands.id:
            for d2 in demands.id:
                if (d1!=d2) and (d1 in demands_fixed) and \
                   (d2 in demands_fixed):
                    Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                         name='Delta_{}_{}'.format(d1, d2), 
                         lb=Deltax[d1, d2], ub=Deltax[d1, d2])
                elif d1!=d2:
                    Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                         name='Delta_{}_{}'.format(d1, d2))
                elif d1==d2:
                    Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                         name='Delta_{}_{}'.format(d1, d2), lb=0, ub=0)
                if mipstart:
                    try:
                        Delta[d1, d2].start = Deltax0[d1, d2]
                    except:
                        pass
#                    
                    
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
            NNN[n] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=10, 
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
            
        # bound for objective
#        model.addConstr(c+Total<=ObjVal+0.5, name='objBound')
        
        # objective
        model.setObjective(c+Total, GRB.MINIMIZE)
            
        # set gurobi parameters
        if len(kwargs):
            for key, value in kwargs.items():
                try:
                    setattr(model.params, key, value)
                except:
                    pass
        model.update()
                
        toc = time.clock()
        self.model_time = toc-tic
                
        model.optimize()
        
        while model.SolCount<1 and num_solve>0:
            if len(kwargs):
                for key, value in kwargs.items():
                    try:
                        setattr(model.params, key, value)
                    except:
                        pass
            model.update()
            model.optimize()
            num_solve -= 1
        
        toc2 = time.clock()
        self.solve_time = toc2-toc
        
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
                        
        # construct solutions for UsageL and Delta
        UsageLx = {}
        for l in self.links.id:
            for d in demands.id:
                if UsageL[l, d].x<0.5:
                    UsageLx[l, d] = 0
                else:
                    UsageLx[l, d] = 1
                       
        Deltax = {}
        for d1 in demands.id:
            for d2 in demands.id:
                if Delta[d1, d2].x <0.5:
                    Deltax[d1, d2] = 0
                else:
                    Deltax[d1, d2] = 1
                          
        Fstartx = {}
        for d in demands.id:
            Fstartx[d] = np.round(Fstart[d].x*2)/2
                   
        Ux = {} # U[a, b] = UsageL[a,b]*Ynode[a,b]
        for l in self.links.id:
            for d in demands.id:
                Ux[l, d] = U[l, d].x
                  
        Irex = {} # 
        for n in self.nodes:
            for d in demands.id:
                Irex[n, d] = Ire[n, d].x
                    
        IIIx = {} # 
        for n in self.nodes:
            for d in demands.id:
                IIIx[n, d] = III[n, d].x
                    
        Ix = {}
        for n in self.nodes:
            Ix[n] = I[n].x
              
        NNNx = {}
        for n in self.nodes:
            NNNx[n] = NNN[n].x
                
        Xx = {}
        for l in self.links.id:
            for d in demands.id:
                Xx[l, d] = X[l, d].x
                  
        Yx = {}
        for n in self.nodes:
            for d in demands.id:
                Yx[n, d] = Y[n, d].x
                      
        Totalx = Total.x
        
        cx = c.x
        
        GNix = {}
        for l in self.links.id:
            for d in demands.id:
                GNix[l, d] = GNi[l, d].x
                    
        G1x={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                G1x[l,d] = G1[l,d].x

        #dvar float GASEws[Links] in 0..10000;
        GASEwsx={}
        for l in self.links.id:
            GASEwsx[l] = GASEws[l].x
             
        #dvar float GNliws[Links][Demands]in 0..10000;
        GNliwsx={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                GNliwsx[l,d] = GNliws[l,d].x
                 
        #dvar float A1[Links][Demands] in 0..10000;
        A1x={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                A1x[l,d] = A1[l,d].x
                 
        #var int UsageL1[Links][Demands][Demands]in 0..1;
        UsageL1x={}# intermedia variables for product
        for l in self.links.id:
            for d1 in demands.id:
                for d2 in demands.id:
                    UsageL1x[l,d1,d2] = UsageL1[l,d1,d2].x
                     
        #dvar float Asenli[Links][Demands]in 0..10000;
        Asenlix={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                Asenlix[l,d] = Asenli[l,d].x
                 
        UseAsenlix = {}
        for l in self.links.id:
            for d in demands.id:
                UseAsenlix[l, d] = UseAsenli[l, d].x
                          
        G1x={}# intermedia variables for product
        for l in self.links.id:
            for d in demands.id:
                G1x[l,d]=G1[l, d].x
        
        
        solutions = {}
        solutions['UsageL'] = UsageLx
        solutions['Fstart'] = Fstartx
        solutions['Delta'] = Deltax
        solutions['U'] = Ux
        solutions['Ire'] = Irex
        solutions['III'] = IIIx
        solutions['I'] = Ix
        solutions['NNN'] = NNNx
        solutions['X'] = Xx
        solutions['Y'] = Yx
        solutions['Total'] = Totalx
        solutions['c'] = cx
        solutions['GASEws'] = GASEwsx
        solutions['GNliws'] = GNliwsx
        solutions['A1'] = A1x
        solutions['UsageL1'] = UsageL1x
        solutions['Asenli'] = Asenlix
        solutions['UseAsenli'] = UseAsenlix
        solutions['G1'] = G1x
        
        return model, solutions, UsageLx, Deltax
    
    def solve_partial_tr(self, demands, previous_solutions, mipstart=False,
                      num_solve=2, FeasibilityTol=1e-9, IntFeasTol=1e-9, 
                      OptimalityTol=1e-9, **kwargs):
        '''Formulate and solve iteratively
        previous_solutions is dict, contains:
            - UsageL from the previous solve, dict
            - Delta from the previous solve, dict
            - demands_added new demands that will be allocated this time, list
            - demands_fixed old demands allocated previously, list, 
                can be empty
            we should make sure that UsageL and Delta contain solutions for 
            all the demands in demands_fixed
        '''
        # process input data
        demands_added = previous_solutions['demands_added'] # new demands
        demands_fixed = previous_solutions['demands_fixed'] # old demands
        UsageLx0 = previous_solutions['UsageL0']
        Deltax0 = previous_solutions['Delta0']
        Fstartx0 = previous_solutions['Fstart0']
        UsageLx = previous_solutions['UsageL']
        Deltax = previous_solutions['Delta']
        ObjVal = previous_solutions['ObjVal']
        demands_all = list(set(demands_added).union(set(demands_fixed)))
#        print(demands_added)
#        print(demands_fixed)
#        print(demands_all)
        n_demands = len(demands_all)
        demands = demands.loc[demands.id.isin(demands_all), :]
        

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
                          
        
        tic = time.clock()
        model = Model('TR')
        
        # define variables
        UsageL = {} # if demand d uses link l
        for l in self.links.id:
            for d in demands.id:
                if d in demands_fixed:
                    UsageL[l, d] = model.addVar(vtype=GRB.BINARY, 
                        name='UsageL_{}_{}'.format(l, d), 
                        lb=UsageLx[l, d], ub=UsageLx[l, d])
                elif d in demands_added:
                    UsageL[l, d] = model.addVar(vtype=GRB.BINARY, 
                        name='UsageL_{}_{}'.format(l, d))
                if mipstart:
                    try:
                        UsageL[l, d].start = UsageLx0[l, d]
                    except:
                        pass
                
        Fstart = {} # the start frequency of demand d
        for d in demands.id:
            Fstart[d] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=bigM1, 
                  name='Fstart_{}'.format(d))
            if mipstart:
                try:
                    Fstart[d].start = Fstartx0[d]
                except:
                    pass
            
        Delta = {} # order between demands
        for d1 in demands.id:
            for d2 in demands.id:
                if (d1!=d2) and (d1 in demands_fixed) and \
                   (d2 in demands_fixed):
                    Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                         name='Delta_{}_{}'.format(d1, d2), 
                         lb=Deltax[d1, d2], ub=Deltax[d1, d2])
                elif d1!=d2:
                    Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                         name='Delta_{}_{}'.format(d1, d2))
                elif d1==d2:
                    Delta[d1, d2] = model.addVar(vtype=GRB.BINARY, 
                         name='Delta_{}_{}'.format(d1, d2), ub=0, lb=0)
                if mipstart:
                    try:
                        Delta[d1, d2].start = Deltax0[d1, d2]
                    except:
                        pass
                    
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
            NNN[n] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=10, 
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
                elif d1==d2:
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
            
        # bound for objective
#        model.addConstr(c+Total<=ObjVal+0.5, name='objBound')
            
        # objective
        model.setObjective(c+Total, GRB.MINIMIZE)
            
        # set gurobi parameters
        if len(kwargs):
            for key, value in kwargs.items():
                try:
                    setattr(model.params, key, value)
                except:
                    pass
        model.update()
                
        toc = time.clock()
        self.model_time = toc-tic
                
        model.optimize()
        
        while model.SolCount<1 and num_solve>0:
            if len(kwargs):
                for key, value in kwargs.items():
                    try:
                        setattr(model.params, key, value)
                    except:
                        pass
            model.update()
            model.optimize()
            num_solve -= 1
        
        toc2 = time.clock()
        self.solve_time = toc2-toc
        
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
                        
        # construct solutions for UsageL and Delta
        UsageLx = {}
        for l in self.links.id:
            for d in demands.id:
                if UsageL[l, d].x<0.5:
                    UsageLx[l, d] = 0
                else:
                    UsageLx[l, d] = 1
                       
        Deltax = {}
        for d1 in demands.id:
            for d2 in demands.id:
                if Delta[d1, d2].x <0.5:
                    Deltax[d1, d2] = 0
                else:
                    Deltax[d1, d2] = 1
                          
        Fstartx = {}
        for d in demands.id:
            Fstartx[d] = np.round(Fstart[d].x*2)/2
                   
        Ux = {} # U[a, b] = UsageL[a,b]*Ynode[a,b]
        for l in self.links.id:
            for d in demands.id:
                Ux[l, d] = U[l, d].x
                  
        Irex = {} # 
        for n in self.nodes:
            for d in demands.id:
                Irex[n, d] = Ire[n, d].x
                    
        IIIx = {} # 
        for n in self.nodes:
            for d in demands.id:
                IIIx[n, d] = III[n, d].x
                    
        Ix = {}
        for n in self.nodes:
            Ix[n] = I[n].x
              
        NNNx = {}
        for n in self.nodes:
            NNNx[n] = NNN[n].x
                
        Xx = {}
        for l in self.links.id:
            for d in demands.id:
                Xx[l, d] = X[l, d].x
                  
        Ynodex = {}
        for n in self.nodes:
            for d in demands.id:
                Ynodex[n, d] = Ynode[n, d].x
                      
        Totalx = Total.x
        
        cx = c.x
        
        solutions = {}
        solutions['UsageL'] = UsageLx
        solutions['Fstart'] = Fstartx
        solutions['Delta'] = Deltax
        solutions['U'] = Ux
        solutions['Ire'] = Irex
        solutions['III'] = IIIx
        solutions['I'] = Ix
        solutions['NNN'] = NNNx
        solutions['X'] = Xx
        solutions['Ynode'] = Ynodex
        solutions['Total'] = Totalx
        solutions['c'] = cx
        
        return model, solutions, UsageLx, Deltax
    
    def scheduler(self, idx, demands, iteration_history, shuffle=False):
        '''Generate demands_fixed and demands_added according to history
            idx is the index of the next step, idx 
            iteration_history contains:
                - step id
                - demands_fixed
                - demands_added
                - solutions
                - UsageLx
                - Deltax
            return:
                - demands_added
                - demands_fixed
                - no_demands, whether there is no demands left and we should 
                    stop
        '''
        
        # calculate timelimit
        timelimit = max(timelimit_baseline, timelimit0*2**(np.floor(idx/n_iter_per_stage)/time_factor) )

        # the first iteration
        if idx==0:
            demands_id = demands.id.as_matrix()
            if shuffle:
                np.random.shuffle(demands_id)
            demands_added = list(demands_id[:n_demands_initial])
            demands_fixed = []
            no_demands = False
            stage_start = False

            return demands_added, demands_fixed, no_demands, timelimit, stage_start

        # how many iterations this set of demands has been solved 
        num_iter_solved = [len(iteration_history[i]['demands_solved']) 
            for i in range(len(iteration_history))]
        num_iter_solved = np.bincount(num_iter_solved)
        # mip gaps
        mipgaps = [iteration_history[i]['model'].mipgap 
            for i in range(len(iteration_history))]

        if num_iter_solved[-1]==n_iter_per_stage:# or mipgaps[-1]<th_mipgap:
            # a set of demands have been solved 10 times
            # or the mip gap is close to zero
            # then finish this round
            
            # check the left demands 
            demands_left = np.setdiff1d(demands.id.as_matrix(),
                        np.array(iteration_history[idx-1]['demands_solved']))
            num_demands_left = len(demands_left)
            if num_demands_left==0:
                # all the demands are solved
                no_demands = True # we don't need to solve anymore
                demands_added = []
                demands_fixed = list(iteration_history[idx-1]['demands_solved'])
                stage_start = False
                
                return demands_added, demands_fixed, no_demands, timelimit, stage_start

            else:
                # there are still some demands left, check the number of left 
                # demands
                num_demands_added = min(n_demands_increment, num_demands_left)
                demands_added = list(demands_left[:num_demands_added])
                demands_fixed = list(iteration_history[idx-1]['demands_solved'])
                no_demands = False
                stage_start = True

                return demands_added, demands_fixed, no_demands, timelimit, stage_start

        else:
            # we are in the middle of a stage, not the first one,
            # so holdout n_demands_holdout demands, add them into the network
            demands_solved = \
                copy.copy(iteration_history[idx-1]['demands_solved'])
            np.random.shuffle(demands_solved)
            n_demands_holdout = int(len(demands_solved)/2)
            demands_added = list(demands_solved[:n_demands_holdout])
            demands_fixed = list(demands_solved[n_demands_holdout:])
            no_demands = False
            stage_start = False

            return demands_added, demands_fixed, no_demands, timelimit, stage_start

    def iterate(self, demands, random_state=0, shuffle=False, mipstart=False, **kwargs):
        '''Solve TR and GN simultaneously, the best solution from TR and GN
            is input to the next iteration of solving (either TR or GN) 
            as the start point
        '''
        tic = time.clock()
        np.random.seed(random_state)
        
        idx = 0
        iteration_history_tr = {}
        iteration_history_tr[idx] = {}
        iteration_history_tr[idx]['step_id'] = idx
        iteration_history_tr[idx]['demands_fixed'] = None
        iteration_history_tr[idx]['demands_added'] = None
        iteration_history_tr[idx]['demands_solved'] = None
        iteration_history_tr[idx]['solutions'] = None
        iteration_history_tr[idx]['UsageLx'] = None
        iteration_history_tr[idx]['Deltax'] = None

        iteration_history_gn = {}
        iteration_history_gn[idx] = {}
        iteration_history_gn[idx]['step_id'] = idx
        iteration_history_gn[idx]['demands_fixed'] = None
        iteration_history_gn[idx]['demands_added'] = None
        iteration_history_gn[idx]['demands_solved'] = None
        iteration_history_gn[idx]['solutions'] = None
        iteration_history_gn[idx]['UsageLx'] = None
        iteration_history_gn[idx]['Deltax'] = None

        # solve the problem for the first time
        # GN model
        demands_added, demands_fixed, _, timelimit, _ = \
            self.scheduler(idx, demands, iteration_history_gn, shuffle)
        demands_tmp = demands.loc[demands.id.isin(demands_added), :]
        model_gn, solutions_gn, UsageLx_gn, Deltax_gn = \
            self.solve_all_gn(demands_tmp, timelimit=timelimit, **kwargs)
        toc_now = time.clock()
        iteration_history_gn[idx]['elapsed_time'] = toc_now-tic
        iteration_history_gn[idx]['demands_fixed'] = demands_fixed
        iteration_history_gn[idx]['demands_added'] = demands_added
        iteration_history_gn[idx]['demands_solved'] = demands_added
        iteration_history_gn[idx]['solutions'] = solutions_gn
        iteration_history_gn[idx]['UsageLx'] = UsageLx_gn
        iteration_history_gn[idx]['Deltax'] = Deltax_gn
        iteration_history_gn[idx]['model'] = model_gn

        # TR model
        model_tr, solutions_tr, UsageLx_tr, Deltax_tr = \
            self.solve_all_tr(demands_tmp, timelimit=timelimit, **kwargs)
        toc_now = time.clock()
        iteration_history_tr[idx]['elapsed_time'] = toc_now-tic
        iteration_history_tr[idx]['demands_fixed'] = demands_fixed
        iteration_history_tr[idx]['demands_added'] = demands_added
        iteration_history_tr[idx]['demands_solved'] = demands_added
        iteration_history_tr[idx]['solutions'] = solutions_tr
        iteration_history_tr[idx]['UsageLx'] = UsageLx_tr
        iteration_history_tr[idx]['Deltax'] = Deltax_tr
        iteration_history_tr[idx]['model'] = model_tr

        idx += 1

        stop_flag = True
        while stop_flag:
            try:
                demands_added, demands_fixed, no_demands, timelimit, stage_start = \
                    self.scheduler(idx, demands, iteration_history_gn, shuffle)

                if no_demands:
                    break

                previous_solutions = {}
                previous_solutions['demands_added'] = demands_added
                previous_solutions['demands_fixed'] = demands_fixed
                # MIPstart
#                if model_tr.ObjVal<model_gn.ObjVal:
#                    previous_solutions['UsageL0'] = UsageLx_tr
#                    previous_solutions['Delta0'] = Deltax_tr
#                    previous_solutions['Fstart0'] = iteration_history_tr[idx-1]['solutions']['Fstart']
#                else:
#                    previous_solutions['UsageL0'] = UsageLx_gn
#                    previous_solutions['Delta0'] = Deltax_gn
#                    previous_solutions['Fstart0'] = iteration_history_gn[idx-1]['solutions']['Fstart']
                
                previous_solutions['UsageL0'] = UsageLx_tr
                previous_solutions['Delta0'] = Deltax_tr
                previous_solutions['Fstart0'] = iteration_history_tr[idx-1]['solutions']['Fstart']
                
                previous_solutions['UsageL'] = UsageLx_tr
                previous_solutions['Delta'] = Deltax_tr
                if stage_start:
                    previous_solutions['ObjVal'] = INF
                else:
                    previous_solutions['ObjVal'] = model_tr.ObjVal

                model_tr, solutions_tr, UsageLx_tr, Deltax_tr = \
                    self.solve_partial_tr(demands, previous_solutions, mipstart=mipstart, 
                                          timelimit=timelimit, **kwargs)
                toc_now = time.clock()
                iteration_history_tr[idx] = {}
                iteration_history_tr[idx]['step_id'] = idx
                iteration_history_tr[idx]['demands_fixed'] = demands_fixed
                iteration_history_tr[idx]['demands_added'] = demands_added
                iteration_history_tr[idx]['demands_solved'] = list(set(demands_fixed).union(set(demands_added)))
                iteration_history_tr[idx]['solutions'] = solutions_tr
                iteration_history_tr[idx]['UsageLx'] = UsageLx_tr
                iteration_history_tr[idx]['Deltax'] = Deltax_tr
                iteration_history_tr[idx]['model'] = model_tr
                iteration_history_tr[idx]['elapsed_time'] = toc_now-tic

#                if model_gn.ObjVal<model_tr.ObjVal:
#                    previous_solutions['UsageL0'] = UsageLx_gn
#                    previous_solutions['Delta0'] = Deltax_gn
#                    previous_solutions['Fstart0'] = iteration_history_gn[idx-1]['solutions']['Fstart']
#                else:
#                    previous_solutions['UsageL0'] = iteration_history_tr[idx]['UsageLx']
#                    previous_solutions['Delta0'] = iteration_history_tr[idx]['Deltax']
#                    previous_solutions['Fstart0'] = iteration_history_tr[idx]['solutions']['Fstart']
                
                previous_solutions['UsageL0'] = UsageLx_gn
                previous_solutions['Delta0'] = Deltax_gn
                previous_solutions['Fstart0'] = iteration_history_gn[idx-1]['solutions']['Fstart']

                previous_solutions['UsageL'] = UsageLx_gn
                previous_solutions['Delta'] = Deltax_gn
                if stage_start:
                    previous_solutions['ObjVal'] = INF
                else:
                    previous_solutions['ObjVal'] = model_gn.ObjVal

                model_gn, solutions_gn, UsageLx_gn, Deltax_gn = \
                    self.solve_partial_gn(demands, previous_solutions, mipstart=mipstart, 
                                          timelimit=timelimit, **kwargs)

                toc_now = time.clock()
                iteration_history_gn[idx] = {}
                iteration_history_gn[idx]['step_id'] = idx
                iteration_history_gn[idx]['demands_fixed'] = demands_fixed
                iteration_history_gn[idx]['demands_added'] = demands_added
                iteration_history_gn[idx]['demands_solved'] = list(set(demands_fixed).union(set(demands_added)))
                iteration_history_gn[idx]['solutions'] = solutions_gn
                iteration_history_gn[idx]['UsageLx'] = UsageLx_gn
                iteration_history_gn[idx]['Deltax'] = Deltax_gn
                iteration_history_gn[idx]['model'] = model_gn
                iteration_history_gn[idx]['elapsed_time'] = toc_now-tic

                idx += 1
                stop_flag = not no_demands

            except:
                try:
                    del previous_solutions
                except:
                    pass
                try:
                    del iteration_history_tr[idx]
                except:
                    pass
                try:
                    del iteration_history_gn[idx]
                except:
                    pass
                gc.collect()

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
    elif hasattr(iteration_history[0]['model'], variable_name):
        var_history = [getattr(iteration_history[i]['model'], variable_name)
            for i in range(len(iteration_history))]

    return var_history

def read_demands(demands_csv, modulation='bpsk'):
    '''read demand csv file'''
    demands = pd.read_csv(demands_csv, header=None)
    demands.reset_index(drop=False, inplace=True)
    demands.columns = ['id', 'source', 'destination', 'data_rates']
    n_demands = demands.shape[0]
    # choose modulation format
    if modulation=='qpsk':
        tr = [(3651-1.25*demands.data_rates[i])/100 for i in range(n_demands)]
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