# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:00:52 2017

@author: misWin

Prepare files for Hebbe, NSF-24, GN vs. TR, bpsk using Xu's algorithm, 
the benchmarks for simulation 3 and 5

Just change the number of iterations per stage to 1 to realize Xu's algorithm
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
num_simulations = 2
n_demands = 20
simulation_name = 'bpsk'
partition = 'economy'
group = 'maite_group'

# resource parameters
cpu_per_job = 1
mem = 64
time_days = 1
time_hours = 0
time_minutes = 0
time_seconds = 0

# input file names
demands_file_template = simulation_name+'_{}.csv'
python_file_template = simulation_name+'_{}.py'
slurm_file_template = simulation_name+'_{}.slurm'
pickle_file_template = simulation_name+'_{}.pkl'

if not os.path.exists(simulation_name):
    os.makedirs(simulation_name)
shutil.copy('milp2.py', simulation_name)
shutil.copy('nsf-24nodes.csv', simulation_name)
shutil.copy('bpsk_TR.csv', simulation_name)
shutil.copy('python_template_bpsk.py', simulation_name)
shutil.copy('batch_template.slurm', simulation_name)
shutil.copy('run_batch_template.py', simulation_name)
os.chdir(simulation_name)

for batch_id in range(num_simulations):
    # generate demands
    network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
    network_cost = network_cost.as_matrix()
    sn = Network(network_cost, modulation='bpsk')
    demands = sn.create_demands(n_demands, modulation='bpsk', low=40, high=100)
    demands_file = demands_file_template.format(batch_id)
    demands.to_csv(demands_file, index=False)
    
    # write python files
    python_src = "python_template_bpsk.py"
    line12 = "batch_id = {}\n".format(batch_id)
    line16 = "demands_file = '{}'\n".format(demands_file)
    pickle_file = pickle_file_template.format(batch_id)
    line24 = "output_file = '{}'\n".format(pickle_file)
    replace_lines = {12:line12, 16:line16, 24:line24}
    python_dst = python_file_template.format(batch_id)
    copy_template(python_src, python_dst, replace_lines)
    
    # write slurm files
    slurm_src = "batch_template.slurm"
    line2 = "#SBATCH --ntasks={}\n".format(cpu_per_job)
    line3 = "#SBATCH --mem={}G\n".format(mem)
    line4 = "#SBATCH --time={}-{}:{}:{}\n".format(time_days, time_hours, 
                        time_minutes, time_seconds)
    line5 = "#SBATCH --job-name={}_{}\n".format(simulation_name, batch_id)
    line6 = "#SBATCH --output={}_{}.stdout\n".format(simulation_name, batch_id)
    line7 = "#SBATCH --error={}_{}.stderr\n".format(simulation_name, batch_id)
    line8 = "#SBATCH --partition={}\n".format(partition)
    line9 = "#SBATCH --account=maite_group\n".format(group)
    line13 = "python {}\n".format(python_dst)
    replace_lines = {2:line2, 3:line3, 4:line4, 5:line5, 6:line6, 7:line7, 
                     8:line8, 9:line9, 13:line13}
    slurm_dst = slurm_file_template.format(batch_id)
    copy_template(slurm_src, slurm_dst, replace_lines)
    
for file in os.listdir(os.curdir):
    change_eol_win2unix(file)
    
os.remove('python_template_bpsk.py')
os.remove('batch_template.slurm')

try:
    os.rename('run_batch_template.py', 'run_batch.py')
except:
    pass
try:
    os.remove('run_batch_template.py')
except:
    pass