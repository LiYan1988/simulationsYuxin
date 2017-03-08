# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 09:49:10 2017

@author: yx4vf

read mat and output csv
"""

import scipy.io as sio
import pandas as pd


for i in range(1, 101):
    file_name = 'demands/demands_14nodes_matlab_{}.mat'.format(i)
    traffic_demands = sio.loadmat(file_name)['m']
    dt = pd.DataFrame({'source':traffic_demands[:, 0]-1, 
    'destination':traffic_demands[:, 1]-1, 'data_rates':traffic_demands[:, 2]})
    dt.reset_index(inplace=True)
    dt = dt[['index', 'source', 'destination', 'data_rates']]
    dt.columns = ['id', 'source', 'destination', 'data_rates']
    bpsk_tr = pd.read_csv('bpsk_TR.csv', header=None)
    bpsk_tr.columns = ['data_rate', 'distance']
    bpsk_tr.distance = bpsk_tr.distance/100
    bpsk_tr.set_index('data_rate', inplace=True)
    tr = [float(bpsk_tr.loc[int(np.round(dt.data_rates[i]))]) 
        for i in dt.id]
    dt['TR'] = tr
    dt.to_csv('demands/demands_14nodes_matlab_{}.csv'.format(i-1))
    