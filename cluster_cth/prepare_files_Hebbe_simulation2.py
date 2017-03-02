# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 19:56:30 2017

@author: misWin

Prepare files for Hebbe, NSF-24, GN vs. TR, qpsk

Note, from bpsk to qpsk, the following should be changed:
    1. when creating demands, indicate the modulation is qpsk, this is for the 
        TR models
    2. when solving milp in python_template_simulation2.py, change the 
        modulation to qpsk, this is for the GN models
    3. since the random seed is fixed, the generated demands should be 
        identical
    4. change template file names
"""

import os
import shutil
import tempfile
from subprocess import call

import pandas as pd
from milp2 import *

np.random.seed(6666)

def copy_template(src, dst, replace_lines):
    destination = open(dst, 'w')
    source = open(src, 'r')
    for l, line in enumerate(source):
        if l in replace_lines.keys():
            destination.write(replace_lines[l])
        else:
            destination.write(line)
    source.close()
    destination.close()
    
def change_eol_win2unix(file_path):
    ftmp, abs_path = tempfile.mkstemp()
    with open(file_path, 'rb') as old_file, open(abs_path, 'wb') as new_file:
        for line in old_file:
            line = line.replace(b'\r\n', b'\n')
            new_file.write(line)
#    
    os.close(ftmp)
    os.remove(file_path)
    os.rename(abs_path, file_path)
    
    
# simulation parameters
num_simulations = 30
n_demands = 50 
simulation_name = 'simulation2_new3' 

# resource parameters
cpu_per_job = 20
time_days = 5
time_hours = 0
time_minutes = 0
time_seconds = 0

# input file names
demands_file_template = simulation_name+'_{}.csv'
python_file_template = simulation_name+'_{}.py'
bash_file_template = simulation_name+'_{}.sh'
pickle_file_template = simulation_name+'_{}.pkl'

if not os.path.exists(simulation_name):
    os.makedirs(simulation_name)
shutil.copy('milp2.py', simulation_name)
shutil.copy('nsf-24nodes.csv', simulation_name)
shutil.copy('bpsk_TR.csv', simulation_name)
shutil.copy('qpsk_TR.csv', simulation_name)
shutil.copy('python_template_simulation2.py', simulation_name)
shutil.copy('batch_template.sh', simulation_name)
shutil.copy('run_batch_template.py', simulation_name)
os.chdir(simulation_name)

for batch_id in range(num_simulations):
    # generate demands
    network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
    network_cost = network_cost.as_matrix()
    sn = Network(network_cost, modulation='qpsk')
    demands = sn.create_demands(n_demands, modulation='qpsk', low=40, high=100)
    demands_file = demands_file_template.format(batch_id)
    demands.to_csv(demands_file, index=False)
    
    # write python files
    python_src = "python_template_simulation2.py"
    line12 = "batch_id = {}\n".format(batch_id)
    line16 = "demands_file = '{}'\n".format(demands_file)
    pickle_file = pickle_file_template.format(batch_id)
    line31 = "output_file = '{}'\n".format(pickle_file)
    replace_lines = {12:line12, 16:line16, 31:line31}
    python_dst = python_file_template.format(batch_id)
    copy_template(python_src, python_dst, replace_lines)
    
    # write bash files
    bash_src = "batch_template.sh"
    line3 = "#SBATCH -J {}_{}\n".format(simulation_name, batch_id)
    line5 = "#SBATCH -n {}\n".format(cpu_per_job)
    line6 = "#SBATCH -t {}-{}:{}:{}\n".format(time_days, time_hours, 
                        time_minutes, time_seconds)
    line7 = "#SBATCH -o {}_{}.stdout\n".format(simulation_name, batch_id)
    line8 = "#SBATCH -e {}_{}.stderr\n".format(simulation_name, batch_id)
    line12 = "pdcp {} $TMPDIR\n".format(python_dst)
    line15 = "pdcp {} $TMPDIR\n".format(demands_file)
    line18 = "python {}\n".format(python_dst)
    replace_lines = {3:line3, 7:line7, 8:line8, 12:line12, 
                     15:line15, 18:line18, 5:line5, 6:line6}
    bash_dst = bash_file_template.format(batch_id)
    copy_template(bash_src, bash_dst, replace_lines)
    
for file in os.listdir(os.curdir):
    change_eol_win2unix(file)
    
os.remove('python_template_simulation2.py')
os.remove('batch_template.sh')

try:
    os.rename('run_batch_template.py', 'run_batch.py')
except:
    pass
try:
    os.remove('run_batch_template.py')
except:
    pass