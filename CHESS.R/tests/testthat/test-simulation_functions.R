test_sets =
  list(
    list(seed = 1),
    
    list(seed = 2),
    
    list(seed = 2, deathrates = 0.2),
    
    list(
      seed = 2,
      birthrates = c(1, 1.2),
      clone_start_times = c(0, 20),
      father = c(0, 0),
      kill_regrow_times = c(0, 0)
    ),
    
    list(
      seed = 2,
      birthrates = c(1, 1.2),
      deathrates = c(0, 0),
      mutation_rates = c(16, 16),
      clone_start_times = c(0, 20),
      father = c(0, 0),
      kill_regrow_times = c(0, 0)
    ),
    
    list(
      seed = 2,
      birthrates = c(1, 1.2),
      deathrates = c(0, 0.1),
      mutation_rates = c(12, 16),
      aggression = c(1, 1),
      push_power = c(1000, 1000),
      clone_start_times = c(0, 20),
      father = c(0, 0),
      kill_regrow_times = c(0, 0)
    ),
    
    list(
      seed = 2,
      birthrates = c(1, 1.2, 1.3),
      deathrates = c(0, 0, 0),
      mutation_rates = c(10, 10, 10),
      aggression = c(1, 1, 1),
      push_power = c(1000, 1000, 1000),
      clone_start_times = c(0, 20, 21),
      father = c(0, 0, 0),
      kill_regrow_times = c(0, 0, 0)
    ),
    
    list(
      seed = 2,
      birthrates = c(1, 1.2, 1.3),
      deathrates = c(0, 0, 0),
      mutation_rates = c(10,10, 10),
      aggression = c(1, 1, 1),
      push_power = c(1000, 1000, 1000),
      clone_start_times = c(0, 20, 21),
      father = c(0, 0, 1),
      kill_regrow_times = c(0, 0, 0)
    )
  )



exp_checksums_samples =
  c(
    "c3aa025b749f4a2db53da7ed63efbffb4a74333e",
    "66f798aaeeacf5ef9ecc24c10b08a806ed8a3739",
    "9b0bdfeb63af6bd487ef65af8a9bb0e498f8a8f4",
    "c06506bb23dafd5b68a2ba8fd10ab5ead52773ab",
    "7c6afd1dca5fbeda9c58d59c11d70e572c4c8924",
    "b77a8e11a0930bcd3512f249812b524f9fcd0200",
    "edabad07c500092e9c13c336f6ffb561965f2c13",
    "d85ddc3ec41ef092cf2ef1bbb55b167b75f08156"
  )

exp_checksums_trees =
  c(
    "05c5bd91faf2e948d542ae39f9abc5a545d7c7c2",
    "16c11af08da6d882cb9c6efc76c11550e517eda9",
    "28bafb88306aaa3b4fae03b943c1ebb34e7f4775",
    "2b573b27e79b5bff37752e9adf12c155b728d2d5",
    "06ed4f3ba58fbe60bde52c380af8d856659934c8",
    "bdb9feeb1f53a297f5346c5b6401ee3f5ed180cf",
    "d0c74608f97a22ccd110ba61d54f2607de7f8ab9",
    "012b13e0807af828912244d79a30a6055761fad9"
  )

tested_sim_samples = lapply(test_sets, do.call, what=SimulateTumourSample)
tested_sim_trees = lapply(test_sets, do.call, what=SimulateTumourTree)



test_that("reproducibility of simulations (trees)", {

  calculate_checksum_trees = 
    function(x) {
      for (i in seq_along(x$tree)) {
        new_labels = gsub("^cell[[:digit:]]+_", "cell_", x$tree[[i]]$tip.label)
        x$tree[[i]]$tip.label = new_labels
        x$tree[[i]]$node.label = NULL
      }
      return(digest::sha1(x))
    }

  checksums_tree = sapply(tested_sim_trees, calculate_checksum_trees)
  expect_equal(as.character(checksums_tree), exp_checksums_trees)
})


test_that("reproducibility of simulations (samples)", {
  
  calculate_checksum_samples = 
    function(x) {
      x[[1]]$mutation_data$id = NULL
      x[[1]]$mutation_data$ccf = NULL
      return(digest::sha1(x))
    }
  
  checksums_samples = sapply(tested_sim_samples, calculate_checksum_samples)
  expect_equal(as.character(checksums_samples), exp_checksums_samples)
})


