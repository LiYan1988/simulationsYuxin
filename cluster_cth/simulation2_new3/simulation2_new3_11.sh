#!/usr/bin/env bash
#SBATCH -A C3SE407-15-3
#SBATCH -p hebbe
#SBATCH -J simulation2_new3_11
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 5-0:0:0
#SBATCH -o simulation2_new3_11.stdout
#SBATCH -e simulation2_new3_11.stderr
module purge 
source ~/.usr_path_grb_py35

pdcp simulation2_new3_11.py $TMPDIR
pdcp milp2.py $TMPDIR
pdcp nsf-24nodes.csv $TMPDIR
pdcp simulation2_new3_11.csv $TMPDIR
cd $TMPDIR

python simulation2_new3_11.py

cp * $SLURM_SUBMIT_DIR

# End script
