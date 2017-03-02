#!/usr/bin/env bash
#SBATCH -A C3SE407-15-3
#SBATCH -p hebbe
#SBATCH -J simulation1_new2_0
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 5-0:0:0
#SBATCH -o simulation1_new2_0.stdout
#SBATCH -e simulation1_new2_0.stderr
module purge 
source ~/.usr_path_grb_py35

pdcp simulation1_new2_0.py $TMPDIR
pdcp milp2.py $TMPDIR
pdcp nsf-24nodes.csv $TMPDIR
pdcp simulation1_new2_0.csv $TMPDIR
cd $TMPDIR

python simulation1_new2_0.py

cp * $SLURM_SUBMIT_DIR

# End script
