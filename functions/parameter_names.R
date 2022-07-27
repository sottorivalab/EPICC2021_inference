# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Style arguments ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

max_sc = 10

parameter_global = c(
  "mutation_rate" = "m", 
  "push_power.1" = "d[push]",
  "deathrates.1" = "mu",
  "region_scale_factor.1" = "d[region]"
)

parameter_order = c(
  "birthrates" = "lambda[%i]",
  "clone_start_times" = "t[%i]",
  "subclone_frac"="N_[%i]/N",
  "per_sample_subclone_frac"="p(SC[%i]|S)"
)

parameter_sc = 
  sprintf(
    rep(parameter_order, each = max_sc),
    rep(seq_len(max_sc), length(parameter_order))
  )


names(parameter_sc) = 
  paste0(
    rep(names(parameter_order), each = max_sc), ".", 
    rep(seq_len(max_sc), length(parameter_order))
  )

params = c(parameter_global, parameter_sc)



model_set_names = c(
  "wgs_only_pushing_to_edge" = "'WGS'",
  "including_lp_pushing_to_edge" = "'WGS & LP'"
)

msi_positiv = c(
  "C536", "C548", "C516", "C518", "C552", "C562"
)

model_string_levels = 
  paste0(
    rep(c("Neutral", "Selection", "Selection x 2", "Selection x 3"), each=2),
    rep(c("", " + Death"), 4)
  )
