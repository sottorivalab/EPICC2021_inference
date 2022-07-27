library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source("setup_environment.R")
add_center_position = CHESS::add_center_position
annotation_from_barcode = THmisc::annotation_from_barcode

options(ignore.interactive=TRUE) # always print progress bar.

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

n_cores = 4
n_top_trees = 0
n_top_sims = 0
force_replot = TRUE
plot_unfished = FALSE
result_source = c("inference_results")
idx = as.numeric(commandArgs(trailingOnly=TRUE)[1])
dump_data = TRUE

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

tree_data = 
  readRDS("input_data/final_tree_set.rds") %>% 
  lapply(function(x){ x$tip.label=gsub(" [(].*", " Added", x$tip.label); return(x)}) %>% 
  lapply(function(x){ x$tip.label=gsub("_D[0-1] Added", "_L1 Added", x$tip.label); return(x)}) %>% 
  lapply(function(x){ x$tip.label=gsub(" Added", "", x$tip.label); return(x)}) # added WGS as LP

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
wh = with(result_dirs_found, (!plots_done | force_replot) & result_set_done & (inference_done | plot_unfished))
result_dirs_found = result_dirs_found$dataset[wh]


for (sim_dir in result_dirs_found) {  
  
  # find number of clonal variants to add to simulations:
  results_this = readRDS(file.path(sim_dir, "result_set.rds"))
  case_id = basename(dirname(sim_dir))
  target_tree = readRDS(file.path(sim_dir, "target.rds"))
  case_tree = ape::keep.tip(tree_data[[case_id]], c("GL", target_tree$tip.label))
  n_clonal = clonal_burden_from_tree(case_tree)
  print(sim_dir)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Create main plots ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
    tryCatch({
      
      n_top = 
        ifelse(
          test = file.exists(file.path(sim_dir, "top_trees.pdf")),
          yes = 0,
          no = n_top_trees
        )
      
      CHESS::plot_smc_results(
        results_this, 
        sim_dir, 
        n_top_sims = n_top_sims, 
        n_top_trees = ,
        n_cores = n_cores,
        labeller_function = labeller_function, 
        plot_tree_function = function(...) plot_tree(lty_lp = 1, alpha_lp = 1, ...), 
        collect_tree_set = 0,
        n_clonal = n_clonal, 
        run_mobster = TRUE, 
        n_muts_mobster = 15000,
        dump_data = dump_data
      )
      
      file.create(file.path(sim_dir, paste0(".replot", .plot_version, ".txt")))
      
    }, error=function(e) {
      print(e)
    })

}

