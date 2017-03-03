#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-0:0:0
#SBATCH --job-name=li_qpsk_0
#SBATCH --output=li_qpsk_0.stdout
#SBATCH --error=li_qpsk_0.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --mem=64G
source ~/.bashrc_grb

python li_qpsk_0.py

