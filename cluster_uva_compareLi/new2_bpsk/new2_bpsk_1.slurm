#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-0:0:0
#SBATCH --job-name=new2_bpsk_1
#SBATCH --output=new2_bpsk_1.stdout
#SBATCH --error=new2_bpsk_1.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --mem=64G
source ~/.bashrc_grb

python new2_bpsk_1.py

