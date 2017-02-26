#!/usr/bin/env bash
#SBATCH -A C3SE407-15-3
#SBATCH -p hebbe
#SBATCH -J test_hebbe_10
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 1-0:0:0
#SBATCH -o test_hebbe_10.stdout
#SBATCH -e test_hebbe_10.stderr
module purge 
source ~/.usr_path_grb_py35

pdcp simulation35xu_10.py $TMPDIR
pdcp milp2_xu.py $TMPDIR
pdcp nsf-24nodes.csv $TMPDIR
pdcp simulation35xu_10.csv $TMPDIR
cd $TMPDIR

python simulation35xu_10.py

cp * $SLURM_SUBMIT_DIR

# End script
