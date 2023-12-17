#!/bin/bash
#SBATCH -M ukko
#SBATCH --job-name=mp_id
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1-20

ROOT_PATH="/home/juhokois/proj/sim_news/local_data/mp_id_"
MANUALLY_ASSIGNED_ID=20231218
NEW_DIR="${ROOT_PATH}${MANUALLY_ASSIGNED_ID}"

mkdir -p $NEW_DIR
module load R/4.2.1

srun Rscript id_sign_slurm.R $NEW_DIR