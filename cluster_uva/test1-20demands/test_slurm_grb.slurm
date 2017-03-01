#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=2-01:00:00
#SBATCH --output=test.stdout
#SBATCH --error=test.stderr
#SBATCH --partition=economy 
#SBATCH --account=optinetly6j

source ~/.bashrc_grb

python test_milp2.py

