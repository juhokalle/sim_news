#!/bin/bash
#SBATCH --job-name=simu
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=1000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-76
#SBATCH -c 1

N_MODS_PER_CORE=8
MANUALLY_ASSIGNED_ID=20220912

srun Rscript b_Rscript4slurm_sim.R $N_MODS_PER_CORE $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID $MANUALLY_ASSIGNED_ID $SLURM_JOB_NAME
