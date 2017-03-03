#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-0:0:0
#SBATCH --job-name=li_qpsk_1
#SBATCH --output=li_qpsk_1.stdout
#SBATCH --error=li_qpsk_1.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --mem=64G
source ~/.bashrc_grb

python li_qpsk_1.py

