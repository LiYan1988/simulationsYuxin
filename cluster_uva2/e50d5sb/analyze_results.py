# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 10:37:28 2017

@author: yx4vf
"""

from milp4 import *
import matplotlib.pyplot as plt

colors = np.array([[0, 0.45, 0.74], [0.85, 0.33, 0.1], [0.93, 0.69, 0.13],
                   [0.49, 0.18, 0.50]])

sn1, iteration_history1 = read_data('e50d5sb_0.pkl')
sn2, iteration_history2 = read_data('e50d4sb_0.pkl')
# spectrum
ctr1 = extract_history(iteration_history1[0], 'c')
cgn1 = extract_history(iteration_history1[1], 'c')
ctr2 = extract_history(iteration_history2[0], 'c')
cgn2 = extract_history(iteration_history2[1], 'c')

# total regeneraters
ttr1 = extract_history(iteration_history1[0], 'Total')
tgn1 = extract_history(iteration_history1[1], 'Total')
ttr2 = extract_history(iteration_history2[0], 'Total')
tgn2 = extract_history(iteration_history2[1], 'Total')

# total circuit
nnntr1 = extract_history(iteration_history1[0], 'NNN')
nnngn1 = extract_history(iteration_history1[1], 'NNN')
nnntr1 = nnntr1[len(nnntr1)-1]
nnngn1 = nnngn1[len(nnngn1)-1]
nnntr1 = [nnntr1[i] for i in range(len(nnntr1))]
nnngn1 = [nnngn1[i] for i in range(len(nnngn1))]

nnntr2 = extract_history(iteration_history2[0], 'NNN')
nnngn2 = extract_history(iteration_history2[1], 'NNN')
nnntr2 = nnntr2[len(nnntr2)-1]
nnngn2 = nnngn2[len(nnngn2)-1]
nnntr2 = [nnntr2[i] for i in range(len(nnntr2))]
nnngn2 = [nnngn2[i] for i in range(len(nnngn2))]


plt.figure(1)
h1, = plt.plot(ctr1, label='TR spectrum new', color=colors[0], linestyle='-.')
h2, = plt.plot(cgn1, label='GN spectrum new', color=colors[1], linestyle='--')
h3, = plt.plot(ctr2, label='TR spectrum old', color=colors[2], linestyle='-.')
h4, = plt.plot(cgn2, label='GN spectrum old', color=colors[3], linestyle='--')
plt.legend(handles=[h1, h2])

plt.figure(2)
h1, = plt.plot(ttr1, label='TR Total new', color=colors[0], linestyle='-.')
h2, = plt.plot(tgn1, label='GN Total new', color=colors[1], linestyle='--')
h3, = plt.plot(ttr2, label='TR Total old', color=colors[2], linestyle='-.')
h4, = plt.plot(tgn2, label='GN Total old', color=colors[3], linestyle='--')
plt.legend(handles=[h1, h2])

plt.figure(4)
ttr1 = np.array(ttr1)
ctr1 = np.array(ctr1)
tgn1 = np.array(tgn1)
cgn1 = np.array(cgn1)

ttr2 = np.array(ttr2)
ctr2 = np.array(ctr2)
tgn2 = np.array(tgn2)
cgn2 = np.array(cgn2)

h1, = plt.plot(np.arange(1, len(ttr1)+1), ttr1+ctr1, label='TR Obj', color=colors[0], linestyle='-.')
h2, = plt.plot(np.arange(1, len(ttr1)+1),tgn1+cgn1, label='GN Obj', color=colors[1], linestyle='--')
h3, = plt.plot(np.arange(1, len(ttr1)+1),ttr1*100, label='TR Total', color=colors[0], linestyle='-.')
h4, = plt.plot(np.arange(1, len(ttr1)+1),tgn1*100, label='GN Total', color=colors[1], linestyle='--')

h5, = plt.plot(np.arange(1, len(ttr2)+1), ttr2+ctr2, label='TR Obj', color=colors[2], linestyle='-.')
h6, = plt.plot(np.arange(1, len(ttr2)+1), tgn2+cgn2, label='GN Obj', color=colors[3], linestyle='--')
h7, = plt.plot(np.arange(1, len(ttr2)+1), ttr2*100, label='TR Total', color=colors[2], linestyle='-.')
h8, = plt.plot(np.arange(1, len(ttr2)+1), tgn2*100, label='GN Total', color=colors[3], linestyle='--')

plt.legend([h1, h2, h3, h4, h5, h6, h7, h8])
#plt.legend([h5, h6, h7, h8])