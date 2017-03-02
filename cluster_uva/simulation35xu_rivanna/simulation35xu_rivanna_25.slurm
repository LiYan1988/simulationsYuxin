#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=1-0:0:0
#SBATCH --job-name=simulation35xu_rivanna_25
#SBATCH --output=simulation35xu_rivanna_25.stdout
#SBATCH --error=simulation35xu_rivanna_25.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python simulation35xu_rivanna_25.py

