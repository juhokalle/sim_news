#!/bin/bash
#SBATCH -M ukko
#SBATCH --job-name=mp_id
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 00:30:00
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-10

ROOT_PATH="/proj/juhokois/sim_news/local_data/mp_id_"
MANUALLY_ASSIGNED_ID=20230908
NEW_DIR="${ROOT_PATH}${MANUALLY_ASSIGNED_ID}"

mkdir -p $NEW_DIR
module load R/4.2.1

srun Rscript b_Rscript4slurm.R $NEW_DIR