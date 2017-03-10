#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5500M
#SBATCH --time=2-0:0:0
#SBATCH --job-name=m4_2
#SBATCH --output=m4_2.stdout
#SBATCH --error=m4_2.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python m4_2.py

