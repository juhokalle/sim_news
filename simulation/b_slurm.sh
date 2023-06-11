#!/bin/bash
#SBATCH -M ukko
#SBATCH --job-name=svarma_simu
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH --ntasks=120
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-96

ROOT_PATH="/proj/juhokois/sim_news/local_data/jobid_"
MANUALLY_ASSIGNED_ID=20230609
NEW_DIR="${ROOT_PATH}${MANUALLY_ASSIGNED_ID}"

mkdir -p $NEW_DIR
module load R/4.2.1

srun Rscript b_Rscript4slurm_sim_new.R $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID $MANUALLY_ASSIGNED_ID $SLURM_ARRAY_TASK_MAX $NEW_DIR