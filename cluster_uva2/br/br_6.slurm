#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=2-0:0:0
#SBATCH --job-name=br_6
#SBATCH --output=br_6.stdout
#SBATCH --error=br_6.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python br_6.py

