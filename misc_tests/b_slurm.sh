#!/bin/bash
#SBATCH -M ukko
#SBATCH --job-name=opt_tests
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH --ntasks=30
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-40

ROOT_PATH="/home/juhokois/proj/sim_news/misc_tests/jobid_"
MANUALLY_ASSIGNED_ID=20231018
NEW_DIR="${ROOT_PATH}${MANUALLY_ASSIGNED_ID}"

mkdir -p $NEW_DIR
module load R/4.2.1

srun Rscript b_Rscript4slurm.R $NEW_DIR