#!/usr/bin/env Rscript
library(dplyr)
library(cowplot)
library(ggplot2)
library(CHESS)
theme_set(theme_cowplot())
options(ignore.interactive=TRUE) # always print progress bar.
source("setup_environment.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

n_cores = 2
force_reload = FALSE
full_chain = TRUE
sim_dirs = c("inference_results")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

sim_dirs = dirname(list.files(sim_dirs, "target[.]rds", recursive = TRUE, full.names = TRUE))
sim_dirs = sim_dirs[file.exists(file.path(sim_dirs, "done.txt"))]

idx = as.numeric(commandArgs(trailingOnly=TRUE)[1])
if (!is.na(idx)) {
  sim_dirs = sim_dirs[idx]
}

for (sim_dir in sim_dirs) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Load data ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  version_file = 
    file.path(
      sim_dir,
      paste0(".results_include_rejected", .result_set_version, ".txt")
    )
  
  results_this =  
    CHESS::update_smc_result_set(
      sim_dir,
      force = force_reload | !file.exists(version_file),
      include_rejected_particles = full_chain,
      n_cores = n_cores
    ) # returns path to rds file or the dataset itself.
  
  file.create(version_file)
  
}
