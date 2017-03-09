#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --time=1-0:0:0
#SBATCH --job-name=bpsk2_0
#SBATCH --output=bpsk2_0.stdout
#SBATCH --error=bpsk2_0.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
source ~/.bashrc_grb

python python_template_simulation35xu.py
python bpsk2_0.py
