#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=2-0:0:0
#SBATCH --job-name=e50d4sb_0
#SBATCH --output=e50d4sb_0.stdout
#SBATCH --error=e50d4sb_0.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python e50d4sb_0.py

