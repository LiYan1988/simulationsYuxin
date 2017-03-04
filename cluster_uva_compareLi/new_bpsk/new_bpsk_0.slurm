#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-0:0:0
#SBATCH --job-name=new_bpsk_0
#SBATCH --output=new_bpsk_0.stdout
#SBATCH --error=new_bpsk_0.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --mem=64G
source ~/.bashrc_grb

python new_bpsk_0.py

