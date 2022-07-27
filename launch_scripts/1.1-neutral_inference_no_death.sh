#!/bin/sh
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=1800MB
#SBATCH --job-name=inf_ntrl
#SBATCH --array=1-32
#SBATCH -o launch_scripts/logs_inference/neutral_modified/inference_neutral_%a_%j.log 

module load R/4.0.5

Rscript ./inference_scripts/01-selection_inference_all_cases_V5.R \
	--N_particles 5000 \
	--M_sim 200 \
	--M_param 1 \
	--idx "$SLURM_ARRAY_TASK_ID" \
	--n_subclones 0 \
	--alpha_param 0.95 \
	--n_cores "$SLURM_NTASKS" \
	--stop_at_n_cells 100000 \
	--edge_distance 0.9 \
	--overdispersion \
	--inference_region_diameter
