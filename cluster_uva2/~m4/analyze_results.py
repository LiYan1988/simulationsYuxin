# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:25:16 2017

@author: yx4vf

analyze results from Rivanna
"""


from milp4 import *
from sys import getsizeof

#sn, iteration_history = read_data('m4_3.pkl')
#
#
##%%
#ctr = extract_history(iteration_history[0], 'c')
#ttr = extract_history(iteration_history[0], 'Total')
#cgn = extract_history(iteration_history[1], 'c')
#tgn = extract_history(iteration_history[1], 'Total')
#
#import matplotlib.pyplot as plt
#plt.figure(1)
#plt.plot(ctr, label='TR')
#plt.plot(cgn, label='GN')
#plt.figure(2)
#plt.plot(ttr, label='TR')
#plt.plot(tgn, label='GN')
#
#nnntr = extract_history(iteration_history[0], 'NNN')
#nnngn = extract_history(iteration_history[1], 'NNN')
#nnntr = [sum(i.values()) for i in nnntr]
#nnngn = [sum(i.values()) for i in nnngn]
#
#plt.figure(3)
#plt.plot(nnntr)
#plt.plot(nnngn)


for i in range(4):
    sn, iteration_history = read_data('m4_{}.pkl'.format(i))
    print('Total time: {}'.format(sn.total_runtime/3600))
    print('Memory: {}'.format(getsizeof(iteration_history)))

