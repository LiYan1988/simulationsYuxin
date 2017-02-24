# -*- coding: utf-8 -*-
"""
Created on Sun Jan 08 09:20:59 2017

@author: celin

simulations for A3
"""

import os
from sdm import *
from shutil import copyfile
import pandas as pd
from milp2 import *

np.random.seed(666)

def copy_template(src, dst, replace_lines):
    destination = open(dst, 'wb')
    source = open(src, 'rU')
    for l, line in enumerate(source):
        if l in replace_lines.keys():
            destination.write(replace_lines[l])
        else:
            destination.write(line)
    source.close()
    destination.close()

num_simulations = 10
n_demands = 5 # for testing 

for batch_id in range(num_simulations):
	demands = sn.create_demands(n_demands, modulation='bpsk', low=40, high=100)
	demands_file = 'nsf24-demand'+batch_id+'.csv'
	demands.to_csv(demands_file)
	src = 'python_template.py'
	line12 = "batch_id = {}\n".format(batch_id)
	line16 = "n_demands = {}\n".format(n_demands)
	replace_lines = {12:line12, 16:line16}
	dst = 'test_hebbe_{}.py'.format(batch_id)
	copy_template(src, dst, replace_lines)