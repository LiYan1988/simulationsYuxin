#!/usr/bin/env bash
#SBATCH -A C3SE407-15-3
#SBATCH -p hebbe
#SBATCH -J batch_template
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH -o batch_template.stdout
#SBATCH -e batch_template.stderr
module purge 
source ~/.usr_path_grb_py35

pdcp python_template.py $TMPDIR
pdcp milp2_xu.py $TMPDIR
pdcp nsf-24nodes.csv $TMPDIR
pdcp demands_template.csv $TMPDIR
cd $TMPDIR

python python_template.py

cp * $SLURM_SUBMIT_DIR

# End script
