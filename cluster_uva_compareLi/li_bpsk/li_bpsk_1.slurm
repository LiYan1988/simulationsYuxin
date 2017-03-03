#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-0:0:0
#SBATCH --job-name=li_bpsk_1
#SBATCH --output=li_bpsk_1.stdout
#SBATCH --error=li_bpsk_1.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --mem=64G
source ~/.bashrc_grb

python li_bpsk_1.py

