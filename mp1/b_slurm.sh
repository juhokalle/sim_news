#!/bin/bash
#SBATCH -M ukko
#SBATCH --job-name=mp_est
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 01:30:00
#SBATCH --ntasks=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-135

ROOT_PATH="/home/juhokois/proj/sim_news/local_data/jobid_"
MANUALLY_ASSIGNED_ID=20231122
NEW_DIR="${ROOT_PATH}${MANUALLY_ASSIGNED_ID}"

mkdir -p $NEW_DIR
module load R/4.2.1

srun Rscript b_Rscript4slurm.R $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID $MANUALLY_ASSIGNED_ID $SLURM_ARRAY_TASK_MAX $NEW_DIR