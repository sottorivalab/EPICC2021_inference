#!/bin/sh
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=5000MB
#SBATCH --job-name=plt_res
#SBATCH -o launch_scripts/logs_other/plot_results_%a_%j.log 

module load R/4.0.5

Rscript ./inference_scripts/06-summary_model_selection.R 
