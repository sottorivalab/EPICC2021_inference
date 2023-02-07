test_that("reproducibility of simulations", {
  
  set.seed(1)
  
  # output dir
  temp_dir = tempfile(pattern = "dir")
  
  # not varied parameters
  params_sim_fixed = list(
    x = 15,
    y = 15,
    z = 1,
    clonal_mutations = 0,
    verbose = FALSE
  )
  
  # particle generator:
  particle_generator = 
    function(x) {
    
    param.matrix =
      rbind(birthrates=1, 
            deathrates=runif(1, 0, 0.49), 
            aggressions=1, 
            push_power=sample(20, 1),
            mutation_rates=10,
            clone_start_times=0,
            kill_regrow_times=0,
            father=-1)
    
    seed = 
      sample.int(1e9, 1)
      
    list(
      param.matrix = param.matrix, 
      seed = seed
    )
    
  }
 
  # Sample generator
  sample_generator = 
    function(sim, seed_samples=NULL) {
      
      n_cells = 8
      
      if (is.null(seed_samples)) seed_samples = sample(.Machine$integer.max, 1)
      set.seed(seed_samples)
      
      # random positions
      # positions = 
      #   with(simulation$params, {
      #     expand.grid(x=seq_len(x), y=seq_len(y), z=seq_len(z)) %>% 
      #       dplyr::mutate(contains_cell = NA)
      #   })

      idx_use = sample(which(sim$coords$type == 1), n_cells)
      
      setup = 
        lapply(idx_use, function(i) {
          center = as.numeric(sim$coords[i,c("x","y","z")])
          diameter = c(1,1,1)
          seed = sample(.Machine$integer.max, 1)
          list(center=center, diameter=diameter, depth=30, purity=1.0, min_vaf=0.05, seed=seed)
        }) %>% magrittr::set_names(. , paste0("S", seq_along(.)))

      return(list(setup=setup, sample_seed=seed_samples))
    }
  
  # target (draw random realization of simulator)
  params_target = particle_generator()
  params_target$param.matrix["push_power",] = 1
  params_target$param.matrix["deathrates",] = 0
  
  target_sim = 
    c(params_sim_fixed, params_target) %>% 
    do.call(what=new_simulation) %>% 
    add_cell_types_CHESS_S3() %>% 
    add_all_samples(., sample_generator(.))
  
  target_tree = 
    with(target_sim, get_equivalent_single_cell_tree(sim, sample_list$setup))
  
  plot(target_sim, sample_tree = FALSE)
  
  upper_bound = lower_bound = params_target$param.matrix
  upper_bound[,1] = c(1, 0.49, 1, 40, 10, 0, 0, -1)
  lower_bound[,1] = c(1, 0, 1, 1, 10, 0, 0, -1)
  
  # SMC test 
  smc_chain = 
    smc_abc_tree(
      const_params = params_sim_fixed, 
      p_gen = particle_generator, 
      s_gen = sample_generator, 
      target_y = target_tree, 
      rho = function(x, y) tree_distance(x, y, "ben"), 
      N = 10, 
      N_T = 8,
      M = 5, 
      M_sim = 5,
      alpha = 0.95,
      epsilon_min = 0, 
      max_steps = 5, 
      n_cores = 1, 
      min_acceptance_rate = 0,
      upper_bound = upper_bound, 
      lower_bound = lower_bound,
      output_dir = file.path(temp_dir, "smc_chain"),
      sample_tree = TRUE, 
      reobserve_resampled = TRUE,
      keep_closest_tree_only = FALSE, 
      allow_increase_in_eps = FALSE, 
      delta_epsilon_min = 0,
      delta_epsilon_wait = 10,
      simulation_manipulator = function(x) { add_cell_types_CHESS_S3(x) }
    )
  
  
  # check produced files
  expected_files = c(
    "smc_chain/options.rds",
    "smc_chain/rejected/S1.rds",
    "smc_chain/rejected/S2.rds",
    "smc_chain/rejected/S3.rds",
    "smc_chain/rejected/S4.rds",
    "smc_chain/rejected/S5.rds",
    "smc_chain/state_step_0.rds",
    "smc_chain/state_step_1.rds",
    "smc_chain/state_step_2.rds",
    "smc_chain/state_step_3.rds",
    "smc_chain/state_step_4.rds",
    "smc_chain/state_step_5.rds"
  )
  
  file_list = list.files(temp_dir, recursive = TRUE)
  expect_equal(file_list, expected_files)
  
 res_data = lapply(file.path(temp_dir, file_list), readRDS)
 expect_known_hash(res_data, "9532b051c7")
 
})


# 
# res = update_smc_result_set(temp_dir)
# if (!"chess_result_set" %in% class(res)) res = readRDS(res)
# 
# plot(cowplot::as_grob(plot_smc_abc_ess(res)))
# plot(cowplot::as_grob(plot_smc_abc_states(res)))
# plot(cowplot::as_grob(posterior_density_plot(res)))
