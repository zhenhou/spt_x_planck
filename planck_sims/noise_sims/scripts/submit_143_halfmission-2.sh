#!/bin/bash
#SBATCH --partition=kicp
#SBATCH --account=kicp
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --job-name=143_mh2
#SBATCH --time=12:00:00

cd /home/zhenhou/Projects/projects/spt_x_planck/planck_sims/noise_sims

mpirun -np 4 ./sxp_sims ini_zone/hfi_143_halfmission-2.ini
