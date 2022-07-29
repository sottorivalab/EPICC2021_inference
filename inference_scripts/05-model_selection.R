library(magrittr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(CHESS)
library(mixR)

options(ignore.interactive=TRUE) # always print progress bar.
theme_set(theme_cowplot())

source("setup_environment.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

n_cores = 18
force_rerun = TRUE
bootstrap_sl_results = TRUE
idx = as.numeric(commandArgs(trailingOnly=TRUE)[1])

N_p = 2000
N_sim = 500
M = 0
N = 2

suffix = "_f2"
particle_file = "result_set.rds"
model_selection_version = "4"


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Find model sets ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# find all model fits
sim_dirs = result_dirs
sim_dirs = sim_dirs[file.exists(file.path(sim_dirs, particle_file))]
names(sim_dirs) = basename(sim_dirs)

# group all model into model sets
dn = dirname(sim_dirs)
model_sets = split(sim_dirs, dn)

if (!is.na(idx)) {
  model_sets = na.omit(model_sets[idx])
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model selection step ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# load model data:
for (i in seq_along(model_sets)) {
  
  # get name of model set:
  model_set_name = names(model_sets)[i]
  model_set = model_sets[[i]]
  print(model_set_name)
  
  marker_files = file.path(model_set, paste0(".model_selection_v", .model_selection_version))

  # rerun model selection
  out_file = file.path(model_set_name, paste0("model_selection_results", suffix, ".rds"))
  if (!file.exists(out_file) | force_rerun | !all(file.exists(marker_files))) {
    
    model_selection = 
      smc_abc_model_selection(
        model_set, 
        suffix = suffix,
        n_cores = n_cores,
        N_p = N_p,
        N = N,
        N_sim = N_sim,
        M = M,
        fname = particle_file,
        title = paste0(basename(model_set_name)),
        keep_trees = FALSE
      )
    
    saveRDS(model_selection, out_file, version = 2)  
    
  } else {
    model_selection = readRDS(out_file)
  }
  

  try({ # get bootstrapped model selection stats
    model_selection_plus_bs = 
      lapply(model_set, function(x) readRDS(file.path(x, paste0("model_selection_data", suffix, ".rds")))) %>% 
      lapply(add_bootstrap_ms_statistics, n=50, n_cores=n_cores) %>% 
      magrittr::set_names(names(model_set))
    out_file = file.path(model_set_name, paste0("model_selection_plus_bs", suffix, ".rds"))
    saveRDS(model_selection_plus_bs, out_file, version = 2)
  })
  

  try({ # get bootstrapped model selection stats with SL
    
    model_selection_plus_bs_syn_lik = 
      lapply(model_set, function(x) readRDS(file.path(x, paste0("model_selection_data", suffix, ".rds")))) %>% 
      lapply(add_bootstrap_ms_statistics_syn, n=50, n_cores=n_cores) %>% 
      magrittr::set_names(names(model_set))
      
    # add delta AIC 
    ids = c("AIC","AIC2")
    for (id in ids) {
      min_ids = do.call(what=rbind, lapply(model_selection_plus_bs_syn_lik, "[[", id))
      min_ids = apply(min_ids, 2, min)
      for (i in seq_along(model_selection_plus_bs_syn_lik)) {
        delta = model_selection_plus_bs_syn_lik[[i]][[id]] - min_ids
        model_selection_plus_bs_syn_lik[[i]][[paste0("delta_", id)]] = delta
      }
    }
    
    ofile = file.path(model_set_name, paste0("model_selection_data_bs", suffix, "_sl.rds"))
    saveRDS(model_selection_plus_bs_syn_lik, ofile)
  })
    

  try({ # plot bootstrapped model selection stats
      
    plot = plot_model_selection_result(
      model_selection_plus_bs_syn_lik,
      title =  paste0("Model selection - ", basename(model_set_name)),
      stats = c("n_param", "neg_log_lik2", "AIC2", "delta_AIC2"),
      alt_labels = c("AIC2"="AIC","neg_log_lik2"="NLL", "delta_AIC2"="Delta~AIC")
    ) + facet_wrap(~L2, nrow=1, scales = "free_y", label = "label_parsed")
      
    plot = plot + 
      geom_hline(data = data.frame(
        L2 = factor("Delta~AIC", levels(plot$data$L2))),
        aes(yintercept = 4),
        linetype = 2
      )

    plot =  plot + # add a polygon showing uncertain range
      geom_polygon(
        data = data.frame(
          x = c(-Inf, Inf, Inf, -Inf, -Inf), y = c(0, 0, 4, 4, 0),
          L2 = factor("Delta~AIC", levels(plot$data$L2)),
          best = TRUE
        ), 
        aes(x = x, y = y), 
        fill="gray90"
      )
    
    plot$layers[c(1,5)] = plot$layers[c(5,1)] # bump added layer down
        
    out_file = file.path(model_set_name, paste0("model_selection_results", suffix, "_syn_lik_bs_alt_aic_only.pdf"))
    ggsave(out_file, plot, height=3.6, width=5.5)
    saveRDS(plot, gsub("[.]pdf$", ".rds", out_file), version=2)
    
  })
    
    
  try({
      
    plot = 
      wrapper_plot_range_of_eps_vs_aic(
        dir = model_set,
        n_cores = n_cores,
        n_points = 50,
        use_sl = TRUE,
        fname = paste0("model_selection_data", suffix, ".rds"),
        bs = 50
      ) + ggtitle(paste0("Model selection - ", basename(model_set_name)))
      
    plot$data =  plot$data %>% 
      dplyr::mutate(L1 = as.character(format_model_string(as.character(L1)))) %>%
      dplyr::filter(variable == "AIC2") %>% 
      dplyr::mutate(variable = "AIC")
      
    plot = plot + 
      ylim(0, NA) + 
      ggtitle(basename(model_set_name)) + 
      scale_x_continuous(n.breaks = 3) + 
      geom_vline(xintercept = model_selection[[1]]$eps_noise, color="gray10", linetype=3) + 
      geom_vline(xintercept = model_selection[[1]]$eps_min_smc, color="gray30", linetype=2)
        
    fname = paste0("model_selection_results", suffix, "_sl_range_of_eps_vs_aic2.pdf")
    out_file = file.path(model_set_name, paste0("model_selection_results", suffix, "_sl_range_of_eps_vs_aic2.pdf"))
    ggsave(out_file, plot, height=2.7, width=4.5)
    saveRDS(plot, gsub("[.]pdf$", ".rds", out_file), version=2)
  })
  
  file.create(marker_files)
}
