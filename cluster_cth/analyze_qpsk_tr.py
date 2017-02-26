# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 20:20:09 2017

@author: misWin

analyze qpsk tr model
"""

import scipy.io as sio
import pandas as pd

qpsk_tr = sio.loadmat('qpsk.mat')['Dis']
idx = list(range(30, 101))
qpsk_tr_ext = [qpsk_tr[0][n-1] for n in idx]
qpsk_TR = pd.DataFrame({'0':idx, '1':qpsk_tr_ext})
qpsk_TR.to_csv('qpsk_TR.csv', header=False, index=False)
qpsk_TR = pd.read_csv('qpsk_TR.csv', header=None)