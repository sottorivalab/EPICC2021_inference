#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1800MB
#SBATCH --job-name=smr_chain
#SBATCH -o launch_scripts/logs_other/summarise_chain_%a_%j.log 

module load R/4.0.5

Rscript ./inference_scripts/02-summary_of_chains.R

