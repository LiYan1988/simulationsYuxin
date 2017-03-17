# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 11:08:44 2017

@author: yx4vf
"""

from milp4 import *
import matplotlib.pyplot as plt


ctr_opt = []
cgn_opt = []
ctr_pro = []
cgn_pro = []

for i in range(20):
    iteration_history_tr, iteration_history_gn = \
        read_data('test_6node_opt{}.pkl'.format(i))
    ctr_opt.append(extract_history(iteration_history_tr, 'c'))
    cgn_opt.append(extract_history(iteration_history_gn, 'c'))
    iteration_history_tr, iteration_history_gn = \
        read_data('test_6node_opt{}.pkl'.format(i))
    ctr_pro.append(extract_history(iteration_history_tr, 'c'))
    cgn_pro.append(extract_history(iteration_history_gn, 'c'))

cgn_opt[19].pop(3)
cgn_opt[18].pop(3)
cgn_pro[18].pop(3)
cgn_pro[19].pop(3)
ctr_opt[19].pop(3)
ctr_opt[18].pop(3)
ctr_pro[18].pop(3)
ctr_pro[19].pop(3)

cgn_pro = np.array(cgn_pro)
cgn_opt = np.array(cgn_opt)
ctr_pro = np.array(ctr_pro)
ctr_opt = np.array(ctr_opt)

plt.figure(1)
plt.plot(cgn_pro.mean(axis=0), label='GN proposed', linestyle='--')
plt.plot(cgn_opt.mean(axis=0), label='GN optimal', linestyle=':')
plt.legend()
plt.show()

plt.figure(2)
plt.plot(ctr_pro.mean(axis=0), label='TR proposed', linestyle='--')
plt.plot(ctr_opt.mean(axis=0), label='TR optimal', linestyle=':')
plt.legend()
plt.show()