#!/bin/bash
#SBATCH --partition=kicp-ht
#SBATCH --account=kicp
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=143_hr2
#SBATCH --time=24:00:00

cd /home/zhenhou/Projects/projects/spt_x_planck/planck_halfring_sims

mpirun -np 1 ./sxp_sims ini_zone/hfi_143_ringhalf_2.ini
