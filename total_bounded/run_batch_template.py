# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 13:34:22 2017

@author: misWin

This is a template file for running batch files on hebbe, python=3.5
"""

import os
from subprocess import call
from shutil import copyfile

for file in os.listdir('.'):
    try:
        extension = file.split('.')[1]
        if extension=='slurm':
            call(['sbatch', file])
    except:
        pass