#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2-01:00:00
#SBATCH --job-name=test_0
#SBATCH --output=test.stdout
#SBATCH --error=test.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --mem=64G
source ~/.bashrc_grb

python python_template_simulation35xu.py

