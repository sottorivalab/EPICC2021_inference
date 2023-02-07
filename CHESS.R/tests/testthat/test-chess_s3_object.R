make_sim_bw_comp = function(x) {
  x$params[c("explore_locally", "return_generation_tree")] = NULL
  
  if ("max_popsize" %in% names(x$params)) 
    if (x$params$max_popsize == 0) 
      x$params$max_popsize = NULL # later added arguments
  
  if ("n_center_estimate" %in% names(x$params)) 
    if (x$params$n_center_estimate == 0) 
      x$params$n_center_estimate = NULL # later added arguments
  
  x$transformation_pos = 
    dplyr::filter(x$transformation_pos, from > -1) # changed behavior
  
  return(x)
}

test_that("constructor of chess S3 object", {
  
  # constructor works for default case
  expect_s3_class(new_simulation(10), "CHESS_S3") # a 10x10x1 simulation, rest default parameters
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # dimension argument
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # wrong parameter for dimensions
  default_size = list(x=10, y=10, z=10)
  expect_s3_class(do.call(what=new_simulation, default_size), "CHESS_S3")
  for (bad_value in list(-1, 0, Inf, -Inf, NA, "A", 0.1, c(5,10), NULL)) {
    for (param in names(default_size)) {
      default_size_mod = default_size
      default_size_mod[param] = list(bad_value)
      expect_error(do.call(what=new_simulation, default_size_mod))
    }
  }
  
  # more bad params for dimensions:
  default_size = list(x=1, y=1, z=1)
  for (bad_value in c(1:2)) {
    for (param in names(default_size)) {
      default_size_mod = default_size
      default_size_mod[param] = list(bad_value)
      expect_error(do.call(what=new_simulation, default_size_mod))
    }
  }
  
  # various good parameter for dimensions
  short_simulation = function(...) {
    suppressWarnings({
      suppressMessages({
        tryCatch({
          R.utils::withTimeout({
            new_simulation(...)
          }, timeout=0.1, onTimeout="silent")
        }, error=function(e) return(NULL))
      }) 
    }) 
    
  warning("1D simulations fail. This must be fixed and tests should be activated.")
  for (x_arg in 4:20) {
    for (y_arg in c(x_arg)) { # for (y_arg in c(1, x_arg)) {
      for (z_arg in c(1, x_arg)) {
        for (i in seq_len(5)) {
            expect_s3_class(short_simulation(x=x_arg, y=y_arg, z=z_arg), "CHESS_S3") 
          }
        }
      }
    }
  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Seed argument
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # check seeding works:
  get_sim_plus_data = function(...) {
    # A little helper function to get the object
    #   including the cell level data, but 
    #   excluding the c++ simulation object
    x = new_simulation(...) %>%
      add_cell_types_CHESS_S3() %>%
      add_cell_mutation_burden_CHESS_S3()
    x$sim = NULL
    return(x)
  }
  
  # test good seed arguments
  sim_d = expect_s3_class(get_sim_plus_data(10, seed=10), "CHESS_S3")
  for (tested_seed in list(1, 12, 123)) {
    tmp = tempfile()
    for (i in seq_len(5)) { # redo simulation 5 times
      sim = get_sim_plus_data(10, seed=tested_seed)
      suppressWarnings(expect_known_output(sim, tmp))
      expect_error(expect_identical(sim, sim_d))
    }
  }
  
  # do we get a random simulation, if we don't set the seed?
  tmp = tempfile()
  suppressWarnings(expect_known_value(get_sim_plus_data(10), tmp)) # create first reference
  for (i in seq_len(10)) {
    expect_error(expect_known_value(get_sim_plus_data(10), tmp, update=TRUE))
  }
  
  # bad arguments
  for (tested_seed in list("A", NA, Inf, -Inf, "1", 0.1, c(1,2))) {
    expect_error(get_sim_plus_data(10, seed=tested_seed))
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Verbosity argument
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # check that switch works
  expect_silent(new_simulation(10, verbose = FALSE))
  expect_output(new_simulation(10, verbose = TRUE))
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # All flag arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  for (arg in c("verbose","explore_locally","record_history")) {
    
    # test some bad values for this argument
    for (arg_value in list(0, 1, "A", NA,"TRUE", 0.1, c(TRUE,FALSE), c(), NULL)) { 
      call_with = magrittr::set_names(list(10, arg_value), c("x", arg))
      expect_error()
    }
    
    # test good values for this argument
    for (arg_value in list(TRUE, FALSE)) { 
      if (arg == "verbose") next() # tested else where and might be noisy ...
      call_with = magrittr::set_names(list(10, arg_value), c("x", arg))
      expect_s3_class(do.call(new_simulation, call_with), "CHESS_S3") # a 10x10x1 simulation, rest default parameters
    }
    
  }
  
  add_data_for_test = function(x) {
    x = add_cell_mutation_burden_CHESS_S3(add_cell_types_CHESS_S3(x))
    x$sim = NULL
    make_sim_bw_comp(x)
  }
  
  expect_known_hash(
    add_data_for_test(
      new_simulation(30, seed = 12)
    ), hash = "af943b91dc")
  
  expect_known_hash(
    add_data_for_test(
      new_simulation(30, seed = 123)
    ), hash = "df2cecbec9")
  
  expect_known_hash(
    add_data_for_test(
      new_simulation(30, seed = 12, clonal_mutations = 500)
  ), hash = "407c914325")
  
  expect_known_hash(
    add_data_for_test(
      new_simulation(30, seed = 12, explore_locally = TRUE)
    ), hash = "36bb64d90d")
  
  expect_known_hash(
    add_data_for_test(
      new_simulation(30, seed = 123, explore_locally = TRUE)
    ), hash = "4f2b6e5a0b")
  
  expect_known_hash(
    add_data_for_test(
      new_simulation(30, seed = 12, clonal_mutations = 500, explore_locally = TRUE)
  ), hash = "67679973ec")
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Try interrupting simulations
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  R.utils::withTimeout(
    expect_warning(new_simulation(200, 200), "Simulation was interrupted."),
    timeout = 0.0001,
    onTimeout = "silent"
  )
  
})

test_that("test snaphot feature", {
  
  # create simulation with snapshots
  snapshot_times = c(1,2,3,10)
  sim = expect_s3_class(new_simulation(20, snapshot_times = snapshot_times), "CHESS_S3") # a 20x20x1 simulation, rest default parameters
  sim_no_snap = expect_s3_class(new_simulation(20), "CHESS_S3") # a 20x20x1 simulation, rest default parameters
  
  # snapshot elements 
  expect_true("snapshots" %in% names(sim))
  expect_false("snapshots" %in% names(sim_no_snap))
  expect_named(sim$snapshots, c("data","time","time_exp"))
  
  # length of snapshot elements
  expect_equal(length(snapshot_times), length(sim$snapshots$data))
  expect_equal(length(snapshot_times), length(sim$snapshots$time))
  expect_equal(length(snapshot_times), length(sim$snapshots$time_exp))
  
  # test animation method
  expect_s3_class(animate(sim), "gganim")
  expect_error(animate(sim_no_snap))
})

test_that("print method", {
  expect_output(print(new_simulation(10, snapshot_times=1:5)))
  expect_output(print(new_simulation(10), "CHESS_S3"))
  expect_output(print(new_simulation(10, snapshot_times=1:5), "CHESS_S3"), regexp = "animate")
})

test_that("cell information methods of chess S3 object", {
  
  sim = new_simulation(20) # a 20x20x1 simulation, rest default parameters
  
  sim_plus_coords = add_coords_CHESS_S3(sim)
  expect_true("coords" %in% names(sim_plus_coords))
  expect_named(sim_plus_coords$coords, c("x","y","z"))
  expect_length(sim_plus_coords$coords$x, 20*20)
  
  
  sim_plus_coords = add_coords_CHESS_S3(sim)
  expect_true("coords" %in% names(sim_plus_coords))
  expect_named(sim_plus_coords$coords, c("x","y","z"))
  expect_length(sim_plus_coords$coords$x, 20*20)
  
  sim_pt = add_cell_mutation_burden_CHESS_S3(add_cell_types_CHESS_S3(sim))
  expect_true("coords" %in% names(sim_pt))
  expect_named(sim_pt$coords, c("x","y","z",'type', 'burden'))
  expect_length(sim_pt$coords$type, 20*20)
  
  expect_equal(sum(sim_pt$coords$type == 1), sim_pt$sim$CellCountsPerType[1])
  expect_gte(min(sim_pt$coords$burden), 0)
  expect_gte(min(sim_pt$coords$type), 0)
  expect_lte(min(sim_pt$coords$type), length(sim_pt$sim$CellCountsPerType))
  
})

test_that("mutation information methods work", {
  
  sim = add_all_samples(new_simulation(20, clonal_mutations = 0), sample_list = "random_cells") # a 20x20x1 simulation, rest default parameters
  
  for (i in seq_along(sim$samples)) {
    # find mutations of samples detected
    wh = sim$samples[[i]]$mutation_data$true_vaf > 0.0 
    mids = sim$samples[[i]]$mutation_data$id[wh]
    # sample mutation burden and gt mutation ids
    args = as.list(sim$samples[[i]]$sample_coord[c(1,2,3)])
    mids_gt = do.call(sim$sim$CellMutations, args)
    m_burden = do.call(sim$sim$CellMutationBurden, args)
    # check that the sets are equal
    expect_true(all(mids %in% mids_gt))
    expect_equal(length(mids_gt), m_burden)
  }
})

test_that("history information methods of chess S3 object", {
  
  sim = expect_s3_class(new_simulation(20, snapshot_times = c(1,2,3,10), record_history=TRUE), "CHESS_S3") # a 20x20x1 simulation, rest default parameters
  expect_true(NROW(sim$sim$History)>1)
  
})

test_that("plotting methods of chess S3 object", {
  
  # create simulation with snapshots
  sim = expect_s3_class(new_simulation(20, snapshot_times = c(1,2,3,10)), "CHESS_S3") # a 20x20x1 simulation, rest default parameters
  sim_no_hist = expect_s3_class(new_simulation(20, snapshot_times = c(1,2,3,10), record_history = FALSE), "CHESS_S3") # a 20x20x1 simulation, rest default parameters
  
  # history plots
  expect_s3_class(plot_history_CHESS_S3(sim), "ggplot")
  expect_s3_class(plot_history_CHESS_S3(sim_no_hist), "ggplot")

  # try plotting the simulation space
  expect_s3_class(plot_space_CHESS_S3(sim), "ggplot")
  expect_s3_class(plot_space_CHESS_S3(sim, alpha_burden=FALSE), "ggplot")
  expect_s3_class(plot_space_CHESS_S3(sim, mm_per_dot=1), "ggplot")
  expect_s3_class(plot_space_CHESS_S3(sim, rotate_by=2), "ggplot")
  expect_s3_class(plot_space_CHESS_S3(sim, rotate_by=2, mirror = TRUE), "ggplot")
})

test_that("sampling method works", {
  
  sim = new_simulation(20) %>% add_cell_types_CHESS_S3()
  
  # try to misuse the sampling method one argument at a time
  bad_args = list(
    center = list(c(1), rep(1, 4), c("A","B","C"), c(NA,1,1), rep(1.1,3)),
    diameter = list(c(1), rep(1, 4), c("A","B","C"), c(NA,1,1), rep(1.1,3)),
    depth = list(-1, Inf, NA, "A"),
    depth_model = list(0, 4, "A", NA, Inf),
    min_reads = list(-0.1, 1.1, -1, "A", Inf, NA),
    min_vaf = list(-0.1, 1.1, "A", Inf, NA),
    seed = list(-0.1, 1.1, "A", Inf, NA),
    purity = list(-0.1, 1.1, "A", Inf, NA)
  )
  
  for (i in seq_along(bad_args)) {
    for (j in seq_along(bad_args[[i]])) {
      args = c(list(x=sim), bad_args[[i]][[j]])
      names(args) = c("x", names(bad_args)[i])
      expect_error(do.call(addSample, args))
    }
  }
  
  # try to use the sampling method one argument at a time
  good_args = list(
    diameter = list(c(2,2,2)),
    depth = list(0, 10, 20),
    depth_model = list(1,2,3),
    min_reads = list(0, 10),
    min_vaf = list(0.0, 0.5, 1.0),
    seed = list(1000, -1000, 0),
    purity = list(0.1, 0.9, 1.0)
  )
  
  for (i in seq_along(good_args)) {
    for (j in seq_along(good_args[[i]])) {
      args = c(list(x=sim, center=c(5,5,1)), good_args[[i]][j])
      names(args) = c("x", "center", names(good_args)[i])
      expect_s3_class(do.call(addSample, args), "CHESS_S3")
    }
  }
  
  
  # check that the add sample method also returns unseen sides with min_reads = 0 and min_vaf = 0,
  # this is important for the multi sample method
  sim_plus_sample = addSample(sim, center = c(10,10,1), depth = 0, min_reads = 0, min_vaf = 0) 
  expect_true(any(sim_plus_sample$samples$S1$mutation_data$depth == 0))


  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # check that results are equal if the order of sampling is changed ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  add_multiple_samples =
    function(s, sample_list) {
      # adds several samples, drop the random mutation ids and remove the simulation object
      
      for (i in seq_along(sample_list))
        s = do.call(addSample, c(list(s), sample_list[[i]]))
      
      for (i in seq_along(s$samples))
        s$samples[[i]]$mutation_data$id = NULL
      
      if (!is.null(names(sample_list)))
        names(s$samples) = names(sample_list)
      
      s$sim = NULL
      
      return(s)
    }
  
  
  sim = expect_s3_class(
    new_simulation(20, clonal_mutations = 50, seed = 10) %>% 
      add_cell_types_CHESS_S3(),
    "CHESS_S3"
  )
  

  # get location of 20 random cells
  set.seed(123)
  wh = sample(which(sim$coords$type > 0), 20)
  idx = as.list(as.data.frame(t(sim$coords[wh,c("x","y","z")])))
  
  
  # 20 random cells to samples:
  samples = lapply(idx, function(x) list(center=x, seed=sample(1e9, 1)))
  names(samples) = paste0("S", seq_along(samples))
  
  
  # try adding sample samples multiple times, same hash?
  expect_known_hash(samples, "e553215ed8")
  for (i in 1:10) {
    sim_plus_samples = make_sim_bw_comp(add_multiple_samples(sim, samples))
    expect_known_hash(sim_plus_samples, "df3f99dd1d")
  }
  
  # permutate the order in which samples are added, same hash?
  for (i in 1:10) {
    ord = sample(seq_along(samples))
    sim_plus_samples = make_sim_bw_comp(add_multiple_samples(sim, samples[ord]))
    sim_plus_samples$samples = sim_plus_samples$samples[order(ord)]
    expect_known_hash(sim_plus_samples, "df3f99dd1d")
  }
  
})

test_that("multiple sampling method works", {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Helper functions ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  check_m_data = function(x) {
    # test that included samples are good
    for (i in seq_along(x$samples)) {
      mdata_i = x$samples[[i]]$mutation_data
      
      # check column names
      expected_cn = c("clone","alt","depth","id","ccf","vaf","cn")
      expect_equal(colnames(mdata_i), expected_cn)
      
      # check mutation data
      wh_not_noise = !grepl("^N", mdata_i$id)
      expect_true(isTRUE(all(mdata_i$clone[wh_not_noise] >= 0)))
      expect_true(isTRUE(all(mdata_i$alt >= 0)))
      expect_true(isTRUE(all(mdata_i$alt <= mdata_i$depth)))
      expect_true(isTRUE(all(!is.na(mdata_i$id))))
      expect_true(isTRUE(all(between(mdata_i$ccf, 0, 1))))
      wh_covered = mdata_i$depth > 0 
      expect_true(isTRUE(all(between(mdata_i$vaf[wh_covered], 0, 1))))
      expect_true(isTRUE(all(mdata_i$cn > 0)))
    }
  }
  
  check_mdata_equal = function(sim) {
    # check that each sample contains data for all detected mutations
    m_lists = lapply(sim$samples, function(s) s$mutation_data$id)
    mids_all = unique(unlist(m_lists)) 
    expect_true(all(sapply(m_lists, function(m) all(mids_all %in% m))))
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Generate simulations ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  sim = expect_s3_class(new_simulation(20, clonal_mutations = 50) %>% add_cell_types_CHESS_S3(), "CHESS_S3")
  
  # get location of 20 cells
  wh = sample(which(sim$coords$type > 0), 20)
  idx = as.list(as.data.frame(t(sim$coords[wh,c("x","y","z")])))

  # sample 20 cells:
  args = lapply(idx, function(x) list(center=x, depth=30, purity=runif(1, 0.25, 0.75), min_reads=5))
  names(args) = NULL
  sim_plus_sample = expect_s3_class(add_all_samples(sim, sample_list = args), "CHESS_S3") # a 20x20x1 simulation, rest default parameters
  sim_plus_sample_noisy = expect_s3_class(add_all_samples(sim, sample_list = args, vaf_noise=1e-4, genome_size = 1e3), "CHESS_S3") 
  
  # as above but will all mutations included (also sites with depth == 0):
  args = lapply(idx, function(x) list(center=x, depth=2, purity=runif(1, 0.25, 0.75), min_reads=0, min_vaf=0))
  sim_plus_sample_all_muts = expect_s3_class(add_all_samples(sim, sample_list = args), "CHESS_S3") # a
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Tests ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  sim_lists = list(sim_plus_sample, sim_plus_sample_noisy, sim_plus_sample_all_muts)
  
  for (s in sim_lists) {
    check_mdata_equal(s)
    check_m_data(s)
  }
  
  # check that all mutations included (i.e. also sites with depth == 0):
  for (i in seq_along(sim_plus_sample_all_muts$samples)) {
    # find mutations of samples detected
    mids = sim_plus_sample_all_muts$samples[[i]]$mutation_data$id
    # sample mutation burden and gt mutation ids
    args = as.list(sim_plus_sample_all_muts$samples[[i]]$sample_coord[c(1,2,3)])
    mids_gt = do.call(sim_plus_sample_all_muts$sim$CellMutations, args)
    # check that the sets are equal
    expect_true(all(mids_gt %in% mids))
  }
  
  
  # try getting mutation table
  vaf = vaf_table_from_samples(sim_plus_sample)
  expect_true(all(!is.na(vaf)))
  expect_true(all(between(unlist(vaf), 0, 1)))
  
  # this does not work, but returns NULL
  expect_null(vaf_table_from_samples(sim))
  
})

test_that("tree retrival methods work", {
  
  sim = 
    expect_s3_class({
      new_simulation(50, clonal_mutations = 50, seed=2)
    }, "CHESS_S3")
  

  # add one low depth sample
  sim_plus_samples = sim %>%  add_all_samples("random_cells")
  sim_plus_samples$sample_list = head(sim_plus_samples$sample_list, 5)
  samples = sim_plus_samples$sample_list
  samples_plus_lp = c(samples, list(list(center=c(25, 25, 1), diameter=c(1, 1, 1), depth=2)))
  sim_plus_samples_and_lp = add_all_samples(sim_plus_samples, samples_plus_lp)
    
  
  #
  tree = get_tree(sim_plus_samples_and_lp, sample_tree=TRUE, add_lowdepth=TRUE)$tree
  tree$node.label = NULL
  expect_s3_class(tree, "phylo")
  expect_equal(length(tree$tip.label), 7)
  expect_equal(ape::write.tree(tree), "(S4:114,(((S1:85,(S3:114,S6:0):10):35,S5:144):54,S2:63):17,GL:0):50;")
  
  tree = get_tree(sim_plus_samples_and_lp, sample_tree=TRUE, add_lowdepth=FALSE)$tree
  tree$node.label = NULL
  expect_s3_class(tree, "phylo")
  expect_equal(length(tree$tip.label), 6)
  expect_equal(ape::write.tree(tree), "(S4:114,(((S1:85,S3:124):35,S5:144):54,S2:63):17,GL:0):50;")
  
  tree = get_tree(sim_plus_samples_and_lp, sample_tree=FALSE, add_lowdepth=TRUE, min_confidence=0, optimize_values=FALSE)$tree
  tree$edge.length = round(tree$edge.length)
  expect_s3_class(tree, "phylo")
  expect_equal(length(tree$tip.label), 7)
  expect_equal(ape::write.tree(tree), "(GL:0,(((S5:144,((S3:110,S6_-Added_confidence-_1:0):14,S1:85):35):54,S2:63):17,S4:114):50):0;")

  tree = get_tree(sim_plus_samples_and_lp, sample_tree=FALSE, add_lowdepth=FALSE)$tree
  expect_s3_class(tree, "phylo")
  expect_equal(length(tree$tip.label), 6)
  expect_equal(ape::write.tree(tree), "((((S5:144,(S3:124,S1:85):35):54,S2:63):17,S4:114):25,GL:25);")

})



test_that("simulation of generation trees work correct", {
  
  check_local_ids_gen_tree =
    function(s) {
      
      mut_ids =
        lapply(seq_len(s$params$x), function(x) {
          lapply(seq_len(s$params$y), function(y) {
            s$sim$CellMutations(x - 1, y - 1, 0)
          }) %>% unlist()
        }) %>% unlist()
      
      unique_mut_ids = sort(strtoi(gsub("X.*", "", unique(mut_ids)), base = 16))
      unique_mut_ids = unique_mut_ids - min(unique_mut_ids) + 1
      expect_equal(unique_mut_ids, seq_along(unique_mut_ids))
    }
  
  
  # constructor works for default case
  sim = expect_s3_class(new_simulation(100, return_generation_tree=TRUE, seed = 123), "CHESS_S3")
  sim_local = expect_s3_class(new_simulation(100, return_generation_tree=TRUE, explore_locally = TRUE, seed = 123), "CHESS_S3")
  
  check_local_ids_gen_tree(sim)
  check_local_ids_gen_tree(sim_local)
})

