library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
source("setup_environment.R")
options(ignore.interactive=TRUE) # always print progress bar.

add_center_position = CHESS::add_center_position
annotation_from_barcode = THmisc::annotation_from_barcode

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

n_cores = 4
force_reexport_posterior_particles = TRUE
result_source = c("inference_results")
load_unfished = FALSE
idx = as.numeric(commandArgs(trailingOnly=TRUE)[1])
n_sample_clones_posterior_particles = 500
last_n_states_posterior_particles = Inf

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Main script ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

result_dirs = file.path(result_source, dirname(list.files(result_source, "target.rds", recursive = TRUE)))
result_dirs = result_dirs[grepl("C[0-9]{3}", basename(dirname(result_dirs)))]
result_dirs_found = print_summary_infos_result_sets(result_dirs, return_list = TRUE)

if (!is.na(idx)) {
  result_dirs_found = result_dirs_found[idx,]
}

# subset to those where plots or particle sets are missing
wh = (!result_dirs_found$posterior_particles_done | force_reexport_posterior_particles) & (result_dirs_found$inference_done | load_unfished) & result_dirs_found$result_set_done
result_dirs_found = result_dirs_found$dataset[wh]

for (sim_dir in result_dirs_found) {  

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Load posterior particle set ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  ifile = file.path(sim_dir, "result_set.rds")
  ofile = file.path(sim_dir, paste0("posterior_particles_v", .posterior_particle_version, ".rds"))
  
  if (file.exists(ifile)) {
    mtime_i = file.mtime(ifile)
    if (!isTRUE(mtime_i <= file.mtime(ofile)) | force_reexport_posterior_particles) {
      
      post_set = 
        ifile %>% 
        load_posterior_particles(last_n_states = last_n_states_posterior_particles) %>% 
        add_sampled_information(n_sample_clones = n_sample_clones_posterior_particles, n_cores = n_cores)
      
      saveRDS(post_set, ofile)
      Sys.setFileTime(ofile, mtime_i)
    }
  }
  
}

