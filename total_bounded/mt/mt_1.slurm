#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=2-0:0:0
#SBATCH --job-name=mt_1
#SBATCH --output=mt_1.stdout
#SBATCH --error=mt_1.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python mt_1.py

