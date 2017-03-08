#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1-0:0:0
#SBATCH --job-name=bpsk_1
#SBATCH --output=bpsk_1.stdout
#SBATCH --error=bpsk_1.stderr
#SBATCH --partition=serial
#SBATCH --account=maite_group

source ~/.bashrc_grb

python bpsk_1.py

