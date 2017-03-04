#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-0:0:0
#SBATCH --job-name=new_bpsk_1
#SBATCH --output=new_bpsk_1.stdout
#SBATCH --error=new_bpsk_1.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --mem=64G
source ~/.bashrc_grb

python new_bpsk_1.py

