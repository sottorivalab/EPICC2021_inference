#!/bin/sh
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=5000MB
#SBATCH --job-name=rld_rplt
#SBATCH --array=1-79
#SBATCH -o launch_scripts/logs_other/reload_and_replot/reload_and_replot_%a_%j.log 

module load R/4.0.5

Rscript ./inference_scripts/03.1-update_result_sets.R "$SLURM_ARRAY_TASK_ID"
Rscript ./inference_scripts/03.2-update_pp_set.R "$SLURM_ARRAY_TASK_ID"
Rscript ./inference_scripts/03.3-replot_results_of_inferences.R "$SLURM_ARRAY_TASK_ID"
Rscript ./inference_scripts/03.4-replot_pp_data.R "$SLURM_ARRAY_TASK_ID"
