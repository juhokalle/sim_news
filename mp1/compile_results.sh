#!/bin/bash
#SBATCH -M ukko
#SBATCH --job-name=svarma_simu
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 00:05:00
#SBATCH --mem=4000
#SBATCH -o rro%a.out
#SBATCH -e rre%a.err
#SBATCH --array=1

ROOT_PATH="tt_empex_"
MANUALLY_ASSIGNED_ID=20230912
FILE_NAME="${ROOT_PATH}${MANUALLY_ASSIGNED_ID}"

module load R/4.2.1

srun Rscript compile_results.R $FILE_NAME