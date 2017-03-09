#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=45G
#SBATCH --time=1-0:0:0
#SBATCH --job-name=bpsk_0
#SBATCH --output=bpsk_0.stdout
#SBATCH --error=bpsk_0.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python bpsk_0.py

