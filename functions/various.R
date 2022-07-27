print_summary_infos_result_sets = function(x, ..., return_list=FALSE) {
  
  unfinished_set = !file.exists(file.path(x, "done.txt"))
  cat(sprintf(
    "Unfinished sets %i/%i (%.2f%%)\n",
    sum(unfinished_set),
    length(unfinished_set),
    mean(unfinished_set) * 100
  ))
  
  missing_res_set = !file.exists(file.path(x, "result_set.rds"))
  not_correct_version = !file.exists(file.path(x, paste0(".results_include_rejected", .result_set_version, ".txt"))) 
  missing_res_set = missing_res_set | not_correct_version
  cat(sprintf(
    "Missing result sets in %i/%i (%.2f%%)\n",
    sum(missing_res_set),
    length(missing_res_set),
    mean(missing_res_set) * 100
  ))
  
  missing_plots = !file.exists(file.path(x, paste0(".replot", .plot_version, ".txt")))
  cat(sprintf(
    "Missing plots in %i/%i (%.2f%%)\n",
    sum(missing_plots),
    length(missing_plots),
    mean(missing_plots) * 100
  ))
  
  missing_pp = !file.exists(file.path(x, paste0("posterior_particles_v", .posterior_particle_version, ".rds")))
  cat(sprintf(
    "Missing posterior particle sets in %i/%i (%.2f%%)\n",
    sum(missing_pp),
    length(missing_pp),
    mean(missing_pp) * 100
  ))
  
  missing_ms_data = !file.exists(file.path(x, paste0(".model_selection_v", .model_selection_version)))
  cat(sprintf(
    "Missing model selection data in %i/%i (%.2f%%)\n",
    sum(missing_ms_data),
    length(missing_ms_data),
    mean(missing_ms_data) * 100
  ))

  if (return_list) {
    res = data.frame(
      dataset = x, 
      inference_done = !unfinished_set,
      result_set_done = !missing_res_set,
      plots_done = !missing_plots,
      posterior_particles_done = !missing_pp,
      model_selection_done = !missing_ms_data
    )
  } else {
    wh = unfinished_set | missing_res_set | missing_plots | missing_pp
    res = x[wh]
  }
  
  invisible(res)
}

