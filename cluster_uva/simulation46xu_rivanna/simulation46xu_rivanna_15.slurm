#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=1-0:0:0
#SBATCH --job-name=simulation46xu_rivanna_15
#SBATCH --output=simulation46xu_rivanna_15.stdout
#SBATCH --error=simulation46xu_rivanna_15.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python simulation46xu_rivanna_15.py

