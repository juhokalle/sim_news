#!/bin/bash
#SBATCH -M ukko
#SBATCH --job-name=svarma_simu
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -n 20
#SBATCH --mem=1000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-100
#SBATCH --constraint=amd

MANUALLY_ASSIGNED_ID=20230407

module load R/4.2.1

srun Rscript b_Rscript4slurm_sim.R $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID $MANUALLY_ASSIGNED_ID $SLURM_ARRAY_TASK_MAX
