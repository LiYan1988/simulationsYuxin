#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2-0:0:0
#SBATCH --job-name=bpsk-serial_1
#SBATCH --output=bpsk-serial_1.stdout
#SBATCH --error=bpsk-serial_1.stderr
#SBATCH --partition=serial
#SBATCH --account=maite_group

source ~/.bashrc_grb

python bpsk-serial_1.py

