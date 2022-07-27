.result_dir = "inference_results" # location of main result files

.expected_inference_set =
  c(LP_and_WGS="abc_smc/including_lp_pushing_to_edge/alma_5000x200x1_alpha_0.95_overdispersed_v3_popsize_100000",
    WGS="abc_smc/wgs_only_pushing_to_edge/alma_5000x200x1_alpha_0.95_overdispersed_v3_popsize_100000"
  )

.expected_model_fits = 
  list(
    LP_and_WGS=c(
      "neutral_no_death_push_regionsize_dispersion", 
      "selection_no_death_push_regionsize_dispersion"
    ), 
    WGS=c(
      "neutral_no_death_push_regionsize_dispersion",
      "selection_no_death_push_regionsize_dispersion"
    )
  )

.optional_model_fits = 
  list(
    LP_and_WGS=c(
      "selection_2_no_death_push_regionsize_dispersion"
    ), 
    WGS=c(
      "selection_2_no_death_push_regionsize_dispersion"
    )
  )

.expected_files_per_case =
  c(
    "best_fit.pdf",
    "region_positions.pdf",
    "posterior_density.png",
    "subclone_fractions.pdf",
    "top_sims.pdf",
    "smc_abc_states.pdf"
  )

.expected_files_per_model = c(
  "best_fit.pdf",
  "region_positions.pdf",
  "posterior_density.png",
  "subclone_fractions.pdf",
  "top_sims.pdf",
  "smc_abc_states.pdf",
  paste0(".replot", .plot_version, ".txt"),
  paste0(".model_selection_v", .model_selection_version),
  "done.txt"
)

.no_content_files = c(
  paste0(".replot", .plot_version, ".txt"),
  paste0(".model_selection_v", .model_selection_version),
  "done.txt"
)

.excluded_cases = c("C519", "C527", "C547")
