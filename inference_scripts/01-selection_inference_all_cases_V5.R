#!/usr/bin/env Rscript
library(CHESS) # simulator package
library(THmisc)
library(gridExtra)
library(cowplot)
library(ggplot2)
theme_set(theme_cowplot())
library(dplyr)
#library(MCMCpack)
library(argparse)
options(ignore.interactive=TRUE) # always print progress bar.
source("setup_environment.R")

output_dir = "inference_results"


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# General arguments ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

tree_file = "input_data/final_tree_set.rds"
particle_collection_dir = NULL 

N_particles = 500
M_sim = 100
M_param = 25
max_steps = 40 
frac_nt = 0.75
sample_tree = TRUE
min_acceptance_rate = 0.05

max_sel_coef = 25
min_delta_eps = 0.01
delta_epsilon_wait = 3

d_method = "sEPD"
inference_push = TRUE
inference_mutation_rate = FALSE
reobserve_resampled = TRUE
allow_increase_in_eps = TRUE
push_to_edge = TRUE

excluded_lp_samples = unlist(readr::read_csv("input_data/excluded_lp_samples.csv"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Sampling setup ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

x = 350 # ~100 colonic crypts per square millimeter  Nguyen H et al (2010). PMID 20689513.
y = 350 # ~colon tumour around 3.5 cm in diameter, i.e. 960 mm^2, 96000 crypts?
z = 1   # diameter of crypt ~ 50um plus 50um of stroma?

depth_wgs = 25   # average sequencing depth
depth_model = 2  # binomial?
depth_lp = 0.5   # average lowpass depth

main_region_diameter = c(25,25,1) # diameter of sampling regions

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Command line arguments
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

parser = ArgumentParser()

parser$add_argument("-i", "--idx", type="integer", default=0, help="Index of sample to use [default %(default)s]", metavar="number")
parser$add_argument("-n", "--n_subclones", type="integer", default=0, help="Number of subclones to use [default %(default)s]", metavar="number")
parser$add_argument("-d", "--inference_deathrate", action='store_true', help="Logical indicating if death rate should be estiamted [default %(default)s]")
parser$add_argument("-a", "--alpha_param", type="double", default=0.90, help="ESS reduction rate [default %(default)s]")
parser$add_argument("-c", "--var_param_clone", type="character", default="birthrate", help="Varied parameter [default %(default)s]")
parser$add_argument("-t", "--test", action='store_true', help="Test mode [default %(default)s]")
parser$add_argument("--inference_region_diameter", action='store_true', default=FALSE, help="Estimate sample region diameter? [default %(default)s]")
parser$add_argument("--no_swap_labels", action='store_true', help="Swap sample labels? [default %(default)s]")
parser$add_argument("--exclude_lp", action='store_true', help="Include LP samples? [default %(default)s]")
parser$add_argument("-x", "--n_cores", type="integer", default=1, help="Number of cores [default %(default)s]")
parser$add_argument("-m", "--method", type="character", default="abc_smc", help="ABC method to use. 'abc_smc' or 'rejection_sampling'. [default %(default)s]")
parser$add_argument("-e", "--epsilon", type="integer", default=1e12, help="Distance for method 'rejection_sampling'. [default %(default)s]")
parser$add_argument("-l", "--N_particles", type="integer", default=N_particles, help="Particles to collect for for method 'rejection_sampling'. [default %(default)s]")
parser$add_argument("-b", "--batch_size", type="integer", default=100, help="Particles to collect for for method 'rejection_sampling'. [default %(default)s]")
parser$add_argument("--overdispersion", action="store_true", help="Estimate overdispersion of mutation rate.")
parser$add_argument("--clip_tip_edges", action="store_true", help="Remove the tip edge of glands?")
parser$add_argument("--sample_from_center", action="store_true", help="Remove the tip edge of glands?")
parser$add_argument("-N", "--stop_at_n_cells", type="integer", default=0, help="Number of cells at stop.")
parser$add_argument("-D", "--edge_distance", type="double", default=0.75, help="")
parser$add_argument("-M", "--M_sim", type="integer", default=100, help="")
parser$add_argument("-P", "--M_param", type="integer", default=25, help="")
parser$add_argument("-p", "--edge_distance_prior_strength", type="double", default=1e9, help="")
parser$add_argument("-z", "--angle_prior_strength", type="double", default=1e9, help="")
parser$add_argument("--no_reject_unexpanded", action="store_true", help="Reject simulations with unexpanded clone?")
parser$add_argument("-Z", "--max_retries_expanded", type="integer", default=3, help="Maximum number of rejections for expanding clone.")


args = parser$parse_args()

cat("Arguments:\n")
for (i in seq_along(args)) {
  cat(paste0("  - ", names(args)[i], ": ", args[[i]]), "\n")
  assign(names(args)[i], args[[i]])
}

if (idx == 0) {
  test = TRUE
  idx = 2
  clip_tip_edges = FALSE
  stop_at_n_cells = 100012
  inference_region_diameter = TRUE
  inference_deathrate = FALSE
  sample_from_center = TRUE
  edge_distance_prior_strength = 20
  angle_prior_strength = 50
  edge_distance = 0.9
  overdispersion = TRUE
}

include_lp = !exclude_lp
swap_labels = !no_swap_labels

stopifnot(method %in% c("abc_smc","rejection_sampling"))
stopifnot(var_param_clone %in% c("birthrate","mutationrate"))

dist_function = (function(m, sl) {
  force(m)
  force(sl)
  
  return(
    function(a, b, ...)
      tree_distance(
        a,
        b,
        method = m,
        swap_labels = sl,
        label_group_function = function(x) {
          paste0(
            sapply(strsplit(x, "_"), "[", 3), "-", 
            substr(sapply(strsplit(x, "_"), "[", 5), 1, 1)
          )
        },
        ...
      )
  )
})(d_method, swap_labels)


if (method == "rejection_sampling") {
  M_param = 1
}

if (clip_tip_edges) {
  tree_modification_function = remove_tips_edge_length
} else {
  tree_modification_function = NULL
}

if (stop_at_n_cells) {
  x = round(sqrt(stop_at_n_cells / pi) * 3)
  y = round(sqrt(stop_at_n_cells / pi) * 3)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Output dirs and cluster specific options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 
result_subdir = paste0(
  ifelse(include_lp, "including_lp", "wgs_only"),
  ifelse(push_to_edge, "_pushing_to_edge", "")
)

opath = file.path(output_dir, method, result_subdir)

fname_suffix = paste0(
  N_particles, "x", M_sim, "x", M_param, "_alpha_", alpha_param,
  ifelse(overdispersion, "_overdispersed_v3", "")
)

fname = paste0("alma_", fname_suffix)

result_dir = file.path(opath, fname)

if (test) {
  
  overdispersion = TRUE
  inference_deathrate = TRUE
  result_dir = paste0(result_dir, "___test___")
  N_particles = 50
  M_sim = 10
  M_param = 1
  max_steps = 4
  
  result_dir = 
    file.path(
      output_dir,
      method,
      result_subdir,
      paste0(
        "/tests_",
        N_particles,
        "_alpha_",
        alpha_param,
        ifelse(overdispersion, "_overdispersed_v3", "")
      )
  )
  
  if (method == "rejection_sampling") {
    N_particles = N_particles * 100
    batch_size = ceiling(N_particles / 10)
  }
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


# trees
tree_data = 
  readRDS(tree_file) %>%  
  lapply(function(x){ x$tip.label=gsub(" [(].*", " Added", x$tip.label); return(x)}) %>% 
  lapply(MLLPT::set_lp_tiplength)

# coverage data
coverage_data = readRDS("input_data/average_coverage.rds")
rownames(coverage_data) = coverage_data$sample_barcode


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Get data for case  ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# load the tree data 
case = names(tree_data)[idx]
print(paste0("IDX: ", idx, " (", case, ")"))
target_tree = tree_data[[case]]
target_tree$tip.label = gsub(" Added", "", gsub("_D[0-1] Added", "_L1 Added", target_tree$tip.label)) # added WGS as LP
tip_labels = target_tree$tip.label[target_tree$tip.label != "GL"]

sample_annotation = 
  annotation_from_barcode(tip_labels, TRUE) %>%
  filter(!sample_type %in% c("L","B"))
  
if (!include_lp) sample_annotation = sample_annotation %>% filter(analyte != "L")
sample_annotation = sample_annotation %>% filter(tissue_type == "cancer")
target_tree = ape::keep.tip(target_tree, c(sample_annotation$sample_barcode, "GL"))
n_clonal = clonal_burden_from_tree(target_tree)
target_tree = target_tree %>% MLLPT:::remove_clonal_variants_tree() # need to rerun

if (!is.null(tree_modification_function)) {
  target_tree = tree_modification_function(target_tree)
}


# create output dir
if (var_param_clone == "birthrate") {
  subdir_name = 
    case_when(
      n_subclones == 0 ~ "neutral", 
      n_subclones == 1 ~ "selection", 
      n_subclones > 1 ~ paste0("selection_", n_subclones))
} else if (var_param_clone == "mutationrate"){
  subdir_name = 
    case_when(
      n_subclones == 0 ~ "neutral", 
      n_subclones == 1 ~ "sc_mutationrate", 
      n_subclones > 1 ~ paste0("sc_mutationrate_", n_subclones))
} else {
  stop("")
}

subdir_name = paste0(subdir_name, ifelse(inference_deathrate, "_death", "_no_death"))
subdir_name = paste0(subdir_name, ifelse(inference_push, "_push", "_no_push"))
subdir_name = paste0(subdir_name, ifelse(inference_region_diameter, "_regionsize", ""))
subdir_name = paste0(subdir_name, ifelse(overdispersion, "_dispersion", ""))
subdir_name = paste0(subdir_name, ifelse(no_reject_unexpanded, "_allow_unexpanded", ""))
result_dir = paste0(result_dir, ifelse(clip_tip_edges, "_clip_tips_v2", ""))
result_dir = paste0(result_dir, ifelse(stop_at_n_cells, paste0("_popsize_", stop_at_n_cells), ""))
out_dir_case = file.path(result_dir, case, subdir_name)


if (inference_region_diameter) {
  main_region_diameter = (main_region_diameter - 1) / (max(main_region_diameter, na.rm=TRUE) - 1) + 1e-12
}


# exist if finished previously
if (file.exists(file.path(out_dir_case, "done.txt"))) {
  cat("=> Previously finished.\n")
  quit(save="no")
}


# Continue chain?
if (file.exists(file.path(out_dir_case, "smc_chain", "options.rds"))) {
  
  cat("Continuing chain.\n")
  
  CHESS::continue_chain(
    out_dir_case, 
    particle_collection_dir = particle_collection_dir,
    n_cores = n_cores,
    delta_epsilon_min = min_delta_eps, # allow reducing eps later on ... 
    delta_epsilon_wait = delta_epsilon_wait,
    max_steps = max_steps, 
    rho = dist_function # updated version of distance function ... 
  )
  
} else {
  
  # save target tree
  dir.create(out_dir_case, FALSE, TRUE)
  if (length(target_tree$tip.label) > 2) {
    plot = plot_tree(target_tree)
    ggsave(file.path(out_dir_case, "tree_target.pdf"), plot, width=5, height=5)
  }
  saveRDS(target_tree, file.path(out_dir_case, "target.rds"), version = 2)
  
  
  # if inference impossible mark as finished.
  if (length(target_tree$tip.label) < 4) {
    file.create(file.path(out_dir_case, "done.txt"))
    cat("=> Not enough samples for analysis.\n")
    quit(save="no")
  }
  
  
  # fixed parameters for all simulations
  params_sim_fixed =
    list(
      x = x,
      y = y,
      z = z,
      clonal_mutations = 0,
      explore_locally = push_to_edge,
      max_popsize = stop_at_n_cells, 
      n_center_estimate = ifelse(sample_from_center, 10000, 0)
    )
  
  
  
  # main sampling regions
  sample_generator = 
    get_sample_generator(
      sample_annot = sample_annotation, 
      edge_distance = edge_distance,
      region_diameter = main_region_diameter, 
      coverage_data = coverage_data, 
      depth_lp = 1, 
      depth_wgs = 1, 
      default_purity = 1, 
      sample_width = c(G=1, B=3),
      find_center = sample_from_center, 
      edge_distance_prior_strength = edge_distance_prior_strength, 
      angle_prior_strength = angle_prior_strength
    )
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Prior functions ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  
  prior_sel = NULL
  prior_t = NULL
  
  
  if (!inference_mutation_rate) {
    prior_mu = function() return(NA)
  } else {
    prior_mu = function() runif(1, 0, 500)
  }
  
  if (!inference_push) {
    prior_push = function()
      return(max(c(x, y, z)) / 2)
  } else {
    get_prior_push = 
      function(x = 1, y = 1, z = 1, ...) {
      (function(max_val)
        function(x)
          round(runif(1, 0, max_val)))(ceiling(max(c(x, y, z)) / 2))
    }
    prior_push = do.call(get_prior_push, params_sim_fixed)
  }
  
  if (!inference_deathrate) {
    prior_dr = function() return(0)
  } else {
    prior_dr = function() runif(1, 0, 0.5)
  }
  
  if (n_subclones & var_param_clone == "birthrate") {
    
    prior_sel = (function(max_s) {
      force(max_s)
      return(function()
        runif(1, 1, max_s))
    })(max_sel_coef)
    
    prior_t = (function(max_t) {
      force(max_t)
      return(function(s)
        floor(runif(length(s), 1, max_t)))
    })(pi * (x / 2) ^ 2 * 0.5)
    
  }
  
  if (overdispersion) {
    prior_dispersion = function() abs(rnorm(1, 0, 2))
  } else {
    prior_dispersion = function() return(NULL)
  }
  
  if (var_param_clone == "mutationrate") {
    
    max_mut_increase = 10
    
    prior_mu = (function(max_increase, prior_mutationrate) {
      force(max_increase)
      force(prior_mutationrate)
      return(function() {
        c(1, runif(n_subclones, 1, max_increase)) * prior_mutationrate()
      })
    })(max_mut_increase, prior_mu)
    
    prior_t = function(s) runif(length(s), 1, 1000)
  }
  
  if (inference_region_diameter) {
    prior_sf_region_size = function() return(runif(1, 15, 25))
  } else {
    prior_sf_region_size = NULL
  }
  
  
  if (method == "rejection_sampling" & case == "C536") {
    prior_push = function() runif(1, 0, 25)
    print(case)
  }
  
  # parameter generator
  particle_generator = # prior
    get_generator_particles(
      prior_mu=prior_mu, 
      prior_push=prior_push, 
      prior_deathrate=prior_dr,
      prior_selection=prior_sel, 
      prior_time=prior_t,
      prior_sf_region_size=prior_sf_region_size,
      n_subclone=n_subclones, 
      prior_dispersion=prior_dispersion
    )
  
  # bounds
  upper_bound = particle_generator()$param.matrix
  upper_bound["push_power", ] = max(c(x,y,z) / 2)
  upper_bound["mutation_rates",] = 500
  upper_bound["deathrates",]  = ifelse(!inference_deathrate, 0, 0.49)
  upper_bound["deathrates",]  = ifelse(!inference_deathrate, 0, 0.49)
  upper_bound["birthrates",]  = max_sel_coef
  upper_bound["clone_start_times",] = floor(pi*(x/2)^2*0.5)
  upper_bound["father",] = seq_len(ncol(upper_bound)) - 2
  
  lower_bound = (particle_generator())$param.matrix
  lower_bound[TRUE] = 0
  lower_bound["birthrates",] = 1
  lower_bound["father",1] = -1 
  lower_bound["father",-1] = 0
  
  if (overdispersion) {
    upper_bound["dispersion",]  = 10
    lower_bound["dispersion",]  = 0
  }
  
  if (inference_region_diameter) {
    upper_bound["region_scale_factor",]  = 25
    lower_bound["region_scale_factor",]  = 15
  }


  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Inference ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  
  if (method == "abc_smc") {
    
    # SMC test 
    smc_chain = 
      CHESS::smc_abc_tree(
        const_params = params_sim_fixed, 
        p_gen = particle_generator, 
        s_gen = sample_generator, 
        target_y = target_tree, 
        rho = dist_function, 
        N = N_particles, 
        N_T = ceiling(frac_nt * N_particles),
        M = M_param, 
        M_sim = M_sim,
        alpha = alpha_param,
        epsilon_min = 0, 
        max_steps = max_steps, 
        delta_epsilon_wait = delta_epsilon_wait,
        n_cores = n_cores, 
        min_acceptance_rate = min_acceptance_rate,
        upper_bound = upper_bound, 
        lower_bound = lower_bound,
        output_dir=file.path(out_dir_case, "smc_chain"),
        sample_tree = sample_tree, 
        reobserve_resampled = reobserve_resampled,
        particle_collection_dir = particle_collection_dir, 
        keep_closest_tree_only = FALSE, 
        allow_increase_in_eps=allow_increase_in_eps, 
        delta_epsilon_min = min_delta_eps, 
        tree_manipulator = tree_modification_function,
        min_frac_clones = ifelse(no_reject_unexpanded, -1, 0.01),
        retries_simulation = max_retries_expanded
      )
    
  } else {
    stop("Unknown sampling method.\n")
  }

}


# Mark complete
file.create(
  file.path(out_dir_case, "done.txt")
)


# Load data
tmp =
  update_smc_result_set(
    out_dir_case,
    force = TRUE,
    include_rejected_particles = TRUE,
    n_cores = ceiling(n_cores * 0.8)
  )

file.create(
  file.path(
    out_dir_case, 
    ".results_include_rejected4.txt"
  )
)


results_this = 
  tryCatch(
    readRDS(tmp),
  error = function(e)
    return(tmp)
  )


# Create plots
tryCatch({
  plot_smc_results(
    results_this,
    out_dir_case, 
    n_top_sims = 0, 
    n_top_trees = 0,
    n_cores = n_cores,
    labeller_function = labeller_function, 
    plot_tree_function = function(...) plot_tree_alt(lty_lp = 1, alpha_lp = 1, ...), 
    collect_tree_set = 0,
    n_clonal = n_clonal, 
    run_mobster = TRUE,
    n_muts_mobster = 30000
  )
}, error=function(e) {
  print(e)
})

