#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=2-0:0:0
#SBATCH --job-name=qr_28
#SBATCH --output=qr_28.stdout
#SBATCH --error=qr_28.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group

source ~/.bashrc_grb

python qr_28.py

