library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
source("setup_environment.R")
options(ignore.interactive=TRUE) # always print progress bar.


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

result_source = c("inference_results")
idx = as.numeric(commandArgs(trailingOnly=TRUE)[1])

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
wh = result_dirs_found$inference_done & result_dirs_found$posterior_particles_done
result_dirs_found = result_dirs_found$dataset[wh]


for (sim_dir in result_dirs_found) {  
  
  pp_file = file.path(sim_dir, paste0("posterior_particles_v", .posterior_particle_version, ".rds"))
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plot subclone fracs ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  tryCatch({
    plot = plot_subclone_fracs(pp_file)
    ofile = file.path(sim_dir, "subclone_fractions.pdf")
    ggsave(ofile, plot, height=5, width=5)
  }, error=function(x) {
    print(x)
    # likely neutral, create empty plot
    plot = ggplot(NULL) + theme(axis.line = element_blank())
    ofile = file.path(sim_dir, "subclone_fractions.pdf")
    ggsave(ofile, plot, height=5, width=5)
  })
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plot posterior of region positions #
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  tryCatch({
    plot = plot_region_positions(pp_file)
    ofile = file.path(sim_dir, "region_positions.pdf")
    ggsave(ofile, plot, height=6.5, width=8)
  }, error=function(x) {
    print(x)
  })
  
}

