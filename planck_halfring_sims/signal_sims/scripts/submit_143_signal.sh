#!/bin/bash
#SBATCH --partition=kicp
#SBATCH --account=kicp
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --job-name=143_signal
#SBATCH --time=12:00:00

cd /home/zhenhou/Projects/projects/spt_x_planck/planck_halfring_sims/signal_sims

mpirun -np 4 ./sxp_sims ini_zone/hfi_143_signal.ini
