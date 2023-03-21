#!/bin/bash
#SBATCH --job-name=svarma_simu
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -c 16
#SBATCH --mem=1000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-167

N_MODS_PER_CORE=15
MANUALLY_ASSIGNED_ID=20230227
SRUN_CPUS_PER_TASK=16

module use /appl/modulefiles/all/ #loads all mds
module load R/4.2.1

srun Rscript b_Rscript4slurm_sim.R $N_MODS_PER_CORE $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID $MANUALLY_ASSIGNED_ID $SLURM_JOB_NAME $SRUN_CPUS_PER_TASK
