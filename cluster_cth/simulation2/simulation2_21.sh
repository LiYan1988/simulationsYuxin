#!/usr/bin/env bash
#SBATCH -A C3SE407-15-3
#SBATCH -p hebbe
#SBATCH -J test_hebbe_21
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 2-0:0:0
#SBATCH -o test_hebbe_21.stdout
#SBATCH -e test_hebbe_21.stderr
module purge 
source ~/.usr_path_grb_py35

pdcp simulation2_21.py $TMPDIR
pdcp milp2.py $TMPDIR
pdcp nsf-24nodes.csv $TMPDIR
pdcp simulation2_21.csv $TMPDIR
cd $TMPDIR

python simulation2_21.py

cp * $SLURM_SUBMIT_DIR

# End script
