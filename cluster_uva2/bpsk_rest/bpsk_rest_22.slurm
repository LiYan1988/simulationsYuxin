#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=2-0:0:0
#SBATCH --job-name=bpsk_rest_22
#SBATCH --output=bpsk_rest_22.stdout
#SBATCH --error=bpsk_rest_22.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python bpsk_rest_22.py

