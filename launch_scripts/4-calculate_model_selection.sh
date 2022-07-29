#!/bin/sh
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=1500MB
#SBATCH --job-name=mdl_sel
#SBATCH --array=1-32
#SBATCH -o launch_scripts/logs_other/model_selection/model_selection_%a_%j.log 

module load R/4.0.5

Rscript ./inference_scripts/05-model_selection.R "$SLURM_ARRAY_TASK_ID"
