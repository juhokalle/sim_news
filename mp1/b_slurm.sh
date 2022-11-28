#!/bin/bash
#SBATCH --job-name=mp_svarma0
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 02:00:00
#SBATCH -c 16
#SBATCH --mem=1000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-113

N_MODS_PER_CORE=3
MANUALLY_ASSIGNED_ID=20221118
SRUN_CPUS_PER_TASK=16

module use /appl/modulefiles/all/ #loads all mds
module load R/4.2.1

srun Rscript b_Rscript4slurm.R $N_MODS_PER_CORE $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID $MANUALLY_ASSIGNED_ID $SLURM_JOB_NAME $SRUN_CPUS_PER_TASK
