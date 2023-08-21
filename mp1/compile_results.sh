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

module load R/4.2.1

srun Rscript compile_results.R