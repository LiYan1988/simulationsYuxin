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
from milp4 import *

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

def change_file(file_path, replace_lines):
    ftmp, abs_path = tempfile.mkstemp()
    with open(file_path, 'r') as old_file, open(abs_path, 'w') as new_file:
        for l, line in enumerate(old_file):
            if l in replace_lines.keys():
                new_file.write(replace_lines[l])
            else:
                new_file.write(line)
#
    os.close(ftmp)
    os.remove(file_path)
    os.rename(abs_path, file_path)

# simulation parameters
num_simulations = 30
n_demands = 50
simulation_name = 'e50d30sb2' #'e{}d{}sb'.format(n_demands, num_simulations) #economy 50 demands 30 samples
partition = 'economy'
group = 'maite_group' # or optinetly6j

# input file names
demands_file_template = simulation_name+'_{}.csv'
python_file_template = simulation_name+'_{}.py'
slurm_file_template = simulation_name+'_{}.slurm'
pickle_file_template = simulation_name+'_{}.pkl'
log_file_template = simulation_name+'_{}.log'

# sbatch parameters
ntasks_per_node = 1
cpus_per_task = 4
mem_per_cpu = 2000
time_days = 2
time_hours = 0
time_minutes = 0
time_seconds = 0

# python parameters
kwargs = {'miphint':True, 'mipfocus':1, 'mipgap':0.001, 'presolve':2,
          'logfile':'log_bpsk.log', 'threads':cpus_per_task}

# milp4 parameters
# objective weight
Nmax = 10
epsilon_total = 1
epsilon_nnn = 0/(Nmax*10)
# modelling parameters
bigM1 = 10**4
bigM2 = 10**4
bigM3 = 2*10**4
# scheduler parameters
n_demands_per_stage = 5
n_iter_per_stage = 10 # 10
timelimit_baseline = 150 # 960
timelimit0 = 20 # 60
time_factor = 1.5


if not os.path.exists(simulation_name):
    os.makedirs(simulation_name)
shutil.copy('milp4.py', simulation_name)
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
    seed = np.random.randint(1000000)
    line10 = "seed = {}\n".format(seed)
    kwargs['logfile'] = log_file_template.format(batch_id)
    line12 = "kwargs = {}\n".format(kwargs.__repr__())
    line14 = "batch_id = {}\n".format(batch_id)
    line18 = "demands_file = '{}'\n".format(demands_file)
    pickle_file = pickle_file_template.format(batch_id)
    line24 = "output_file = '{}'\n".format(pickle_file)
    replace_lines = {10:line10, 12:line12, 14:line14, 18:line18, 24:line24}
    python_dst = python_file_template.format(batch_id)
    copy_template(python_src, python_dst, replace_lines)

    # write slurm files
    slurm_src = "batch_template.slurm"
    line2 = "#SBATCH --ntasks-per-node={}\n".format(ntasks_per_node)
    line3 = "#SBATCH --cpus-per-task={}\n".format(cpus_per_task)
    line4 = "#SBATCH --mem-per-cpu={}M\n".format(mem_per_cpu)
    line5 = "#SBATCH --time={}-{}:{}:{}\n".format(time_days, time_hours, time_minutes, time_seconds)
    line6 = "#SBATCH --job-name={}_{}\n".format(simulation_name, batch_id)
    line7 = "#SBATCH --output={}_{}.stdout\n".format(simulation_name, batch_id)
    line8 = "#SBATCH --error={}_{}.stderr\n".format(simulation_name, batch_id)
    line9 = "#SBATCH --partition={}\n".format(partition)
    line10 = "#SBATCH --account=maite_group\n".format(group)
    line14 = "python {}\n".format(python_dst)
    replace_lines = {2:line2, 3:line3, 4:line4, 5:line5, 6:line6, 7:line7,
                     8:line8, 9:line9, 10:line10, 14:line14}
    slurm_dst = slurm_file_template.format(batch_id)
    copy_template(slurm_src, slurm_dst, replace_lines)

# write milp4.py
replace_lines = {}
replace_lines[24] = "Nmax = {}\n".format(Nmax)
replace_lines[25] = "epsilon_total = {}\n".format(epsilon_total)
replace_lines[26] = "epsilon_nnn = {}\n".format(epsilon_nnn)
replace_lines[29] = "bigM1 = {}\n".format(bigM1)
replace_lines[30] = "bigM2 = {}\n".format(bigM2)
replace_lines[31] = "bigM3 = {}\n".format(bigM3)
replace_lines[34] = "n_demands_per_stage = {}\n".format(n_demands_per_stage)
replace_lines[35] = "n_iter_per_stage = {}\n".format(n_iter_per_stage)
replace_lines[36] = "timelimit_baseline = {}\n".format(timelimit_baseline)
replace_lines[37] = "timelimit0 = {}\n".format(timelimit0)
replace_lines[38] = "time_factor = {}\n".format(time_factor)
change_file('milp4.py', replace_lines)

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