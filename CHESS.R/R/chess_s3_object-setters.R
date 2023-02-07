#' Adds coordinate data to a CHESS S3 object.
#'
#' @param x A CHESS_S3 object.
#'
#' @return The CHESS_S3 object with the added data.
#' @export
#'
#'@examples add_coords_CHESS_S3(new_simulation(x=20, y=20)) # 20x20 simulation
add_coords_CHESS_S3 = function(x) {
  
  checkmate::assertClass(x, "CHESS_S3")
  
  if (!"coords" %in% names(x)) {
    x$coords =
      with(x$params, { 
        data.frame(expand.grid(x = seq_len(x), y = seq_len(y), z = seq_len(z)))
      })
  }
  
  return(x)
  
}


#' Adds cell type data to a CHESS S3 object.
#'
#' @param x A CHESS_S3 object.
#'
#' @return The CHESS_S3 object with the added data.
#' @export
#
#'@examples add_cell_types_CHESS_S3(new_simulation(x=20, y=20)) # 20x20 simulation
add_cell_types_CHESS_S3 = function(x) {
  
  checkmate::assertClass(x, "CHESS_S3")
  
  x = add_coords_CHESS_S3(x)
  
  if (!"type" %in% names(x$coords)) {
    get_type = function(d) x$sim$CellType(d[1], d[2], d[3])
    x$coords$type = apply(x$coords-1, 1, get_type)
  }
  
  x$transformation_pos = x$sim$Transformations
  
  return(x)
}


#' Adds each cells mutation burden to a CHESS S3 object.
#'
#' @param x A CHESS_S3 object.
#'
#' @return The CHESS_S3 object with the added data.
#' @export
#'
#'@examples add_cell_mutation_burden_CHESS_S3(new_simulation(x=20, y=20)) # 20x20 simulation
add_cell_mutation_burden_CHESS_S3 = function(x) {
  
  checkmate::assertClass(x, "CHESS_S3")
  
  x = add_coords_CHESS_S3(x)
  
  if (!"burden" %in% names(x$coords)) {
    x$coords$burden = 
      apply(x$coords-1, 1, function(d) x$sim$CellMutationBurden(d[1],d[2],d[3]))
  }
  
  return(x)
}


#' Adds cell type data for a subset of cells to a CHESS S3 object.
#'
#' @param x A CHESS_S3 object.
#'
#' @return The CHESS_S3 object with the added data.
#'
#'@examples add_random_types(new_simulation(x=20, y=20)) # 20x20 simulation
add_random_types = function(x, n=1000) {
  
  # init cell list
  x = add_coords_CHESS_S3(x)
  x$coords$type_rnd = NA
  
  # add random types
  if ("type" %in% colnames(x$coords)) {
    x$coords$type_rnd = x$coords$type
  } else {
    get_type = function(d) x$sim$CellType(d[1], d[2], d[3])
    wh = sample.int(NROW(x$coords), min(c(NROW(x$coords), n)))
    x$coords$type_rnd[wh] =  apply(x$coords[wh,]-1, 1, get_type)
  }
  
  return(x)
}


#' Estimates the center position of the simulations
#'
#' @param x A CHESS_S3 object.
#'
#' @return The CHESS_S3 object with the center added
#' @export
#'
#'@examples add_center_position(new_simulation(x=20, y=20)) # 20x20 simulation
add_center_position = function(x, n=1000) {
  
  if ("estimated_center" %in% names(x)) 
    return(x)
  
  rng_state = .GlobalEnv$.Random.seed
  
  set.seed(x$params$seed)
  x = add_random_types(x, n=n)
  
  x$estimated_center =
    with(
      x$coords[which(x$coords$type_rnd > 0),], 
      data.frame(
        x = mean(x), 
        y = mean(y), 
        z = mean(z)
      )
    )
  
  if (!is.null(rng_state)) {
    .GlobalEnv$.Random.seed = rng_state
  } else {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  
  return(x)
}



#' Adds the simulation history to a CHESS S3 object.
#'
#' @param x A CHESS_S3 object.
#'
#' @return The CHESS_S3 object with the added data.
#' @export
#'
#'@examples add_history_CHESS_S3(new_simulation(x=20, y=20, record_history = TRUE)) # 20x20 simulation
add_history_CHESS_S3 = function(x) {
  
  checkmate::assertClass(x, "CHESS_S3")
  
  if (!"history" %in% names(x)) {
    x$history = x$sim$History
  }
  
  return(x)
}


#' Adds local bulk sample entropy to object
#'
#' @param s 
#' @param radius 
#' @param min_vaf 
#'
#' @return
#' @export
#'
#' @examples add_local_entropy(new_simulation(20), radius=c(2,2,0))
add_local_entropy =
  function(s, radius = c(5, 5, 0), min_vaf=0.05) {
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Check arguments
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    al = checkmate::makeAssertCollection()
    checkmate::assertClass(s, "CHESS_S3")
    checkmate::assertIntegerish(radius, any.missing = FALSE, null.ok = FALSE, len = 3, lower=0, add=al)
    checkmate::assertNumeric(min_vaf, any.missing = FALSE, null.ok = FALSE, lower=0, upper=1, len = 1, add=al)
    checkmate::reportAssertions(al)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Main code
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    cat("Adding local mutation entropy data.\n")
    
    depth = 10000 # depth for the estimation of ccf
    s = add_cell_types_CHESS_S3(s)
    
    # construct a list of arguments for the sample function
    centers = s$coords[, c("x", "y", "z")]
    pos = t(rbind(t(centers) - radius, t(centers) + radius) - 1)
    pos_arg_list = split(pos, seq_len(NROW(pos)))
    const_args = list(1e5, 3, 1, min_vaf, sample.int(1e9, 1), 1.0)
    args = lapply(pos_arg_list, function(x) c(as.list(x), const_args))
    
    # label of the new entropy column
    cn_new = tail(make.unique(c(colnames(s$coords), "entropy")), 1)
    
    # sample tumor and estimate entropy value 
    s$coords[,cn_new] =
      pbapply::pbsapply(args, function(x) {
        sample = do.call(s$sim$TakeSample, x)$mutation_data
        if (NROW(sample) == 0) return(NA)
        ccf = sample$alt / sample$depth * 2; ccf[ccf>1] = 1
        -sum((ccf * log(ccf)) / log(length(ccf)))
      })
    
    return(s)
  }


#' Adds a sample to the simulation object
#'
#' @param x CHESS S3 object.
#' @param center Center of sample position (1-based index).
#' @param diameter Diameter of sample position(1-based index).
#' @param depth Average sequencing depth (default 100x). 
#' @param depth_model Sequencing model.
#' @param min_reads Minimum number of reads to keep mutation (default: 1).
#' @param min_vaf Minimum VAF to keep mutation (default 0.05)..
#' @param seed Seed. Note: A value of 0 indicates to not set a seed for the sample. This will use the state of the random number generator at the time.
#' @param purity Sample purity (default: 1.0, pure sample).
#' @param ... 
#'
#' @return CHESS S3 object with the added sample.
#'
#'@examples add_sample(new_simulation(x=20, y=20, seed=1), center=c(10,10)) # 20x20 simulation
#'
add_sample = function(x, center, diameter=rep(1, length(center)), depth=100, depth_model=1, min_reads=2, min_vaf=0.05, seed=1, purity=1.0, ...) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assertIntegerish(center, any.missing = FALSE, null.ok = FALSE, min.len = 2, max.len = 3, lower=0,  add=al)
  checkmate::assertIntegerish(diameter, any.missing = FALSE, null.ok = FALSE, min.len = 2, max.len = 3, lower=1,  add=al)
  checkmate::assertTRUE(length(diameter) == length(center), add=al)
  checkmate::assertDouble(depth, any.missing = FALSE, null.ok = FALSE, len = 1, lower=0, finite = TRUE, add=al)
  checkmate::assertIntegerish(depth_model, any.missing = FALSE, null.ok = FALSE, len = 1, lower=0, upper = 3, add=al)
  checkmate::assertIntegerish(min_reads, any.missing = FALSE, null.ok = FALSE, len = 1, lower=0, add=al)
  checkmate::assertDouble(min_vaf, any.missing = FALSE, null.ok = FALSE, len = 1, lower=0, upper=1, add=al)
  checkmate::assertIntegerish(seed, any.missing = FALSE, null.ok = FALSE, len = 1, add=al)
  checkmate::assertDouble(purity, any.missing = FALSE, null.ok = FALSE, len = 1, lower=0, upper=1, add=al)
  
  checkmate::reportAssertions(al)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Fix alternative arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  if (length(center) == 2) center = c(center, 1) # add missing z coordinate
  if (length(diameter) == 2) diameter = c(diameter, 1) # add missing z coordinate
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main function code
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  
  if (!"samples" %in% names(x)) {x[["samples"]] = list()}
  
  coords = ceiling(c(center - diameter/2, center + diameter/2 - 1)) - 1 # zero indexed
  
  other_args =
    list(
      depth = depth,
      depth_model = depth_model,
      min_reads = min_reads,
      min_vaf = min_vaf,
      seed = seed,
      purity = purity
    )
  
  
  args = c(sample_coord = as.list(coords), other_args)
  sampled_data = do.call(x$sim$TakeSample, args)
  
  # insert vaf values:
  sampled_data$mutation_data = 
    sampled_data$mutation_data %>% 
    dplyr::mutate(vaf = alt / depth) %>% 
    dplyr::mutate(cn = 2)
  
  
  # append data 
  sample_result = list(c(list(sample_coords=coords), other_args, sampled_data))
  names(sample_result) = paste0("S", length(x$samples)+1)
  x[["samples"]] = c(x[["samples"]], sample_result)
  
  invisible(x)
}


#' Add multiple samples to a CHESS S3 object
#' @param ... 
#'
#' @describeIn add_sample Add multiple samples to a CHESS S3 object
addSample = function(...) {
  add_sample(...) # old name, keep but don't document
}


#' Add multiple samples to a CHESS S3 object. 
#' 
#' Note: Old samples will be removed. Mutations that are detected in a other sample
#' will also be keeped in all other samples. 
#'
#' @param s CHESS S3 object.
#' @param sample_list A list of samples. Each element must contain a list with arguments that can be passed to the add_sample function (see ?add_sample for details). Try 'random_cells' for testing. 
#'
#' @return CHESS S3 object with all the samples added.
#' 
#' @examples add_all_samples(new_simulation(20, seed=1), "random_cells") # a small 20x20x1 simulation
#' 
add_all_samples = function(s, sample_list, vaf_noise=0, genome_size=3.9e9, ...) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Alternative values for optional arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # no sample list passed, use random cells instead.
  if (isTRUE(all.equal("random_cells", sample_list))) {
    s = add_cell_types_CHESS_S3(s)
    set.seed(s$params$seed)
    
    # get location of up to 20 random cells.
    wh_cells = which(s$coords$type > 0)
    wh_sampled = sample(wh_cells, min(c(20, length(wh_cells))))
    pos_list = as.list(as.data.frame(t(s$coords[wh_sampled,c("x","y","z")])))
    sample_list = lapply(pos_list, function(x) list(center=x, diameter=rep(1, length(x))))
    names(sample_list) = NULL
  }
  
  # sample list can come with optional properties,
  # keep the main setup in this case
  sample_list_orig = sample_list
  if ("setup" %in% names(sample_list)) {
    sample_list = sample_list$setup
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  checkmate::assertClass(s, "CHESS_S3", add=al)
  checkmate::assertList(sample_list, null.ok = TRUE, min.len = 0, types="list",  add=al)
  checkmate::assertNumeric(vaf_noise, lower=0, upper=1, any.missing = FALSE, len=1,  null.ok = FALSE, add=al)
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main function code
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  s$samples = NULL # remove all old samples
  
  if (all(c("region_diameter","region_center") %in% names(sample_list_orig))) {
    s$sample_regions = sample_list_orig[c("region_diameter","region_center")]
  }
  s$sample_list = sample_list_orig
  
  if ("setup" %in% names(sample_list_orig)) {
    sample_list = sample_list_orig$setup
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Add all main samples
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # add all samples
  for (i in seq_along(sample_list)) {
    
    # modify arguments
    arg = c(list(x=s), sample_list[[i]])
    arg[c("min_reads","min_vaf")] = list(0, 0.0)
    
    # add sample
    s = do.call(add_sample, arg)
    
    # set name
    if (!is.null(names(sample_list))){
      names(s$samples)[length(s$samples)] = names(sample_list)[i]
    } else {
      names(s$samples)[length(s$samples)] = paste0("S", i)
    }
  }
  names(s$samples) = make.unique(names(s$samples))
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Add pure noise data 
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  n_noisy_sites = 0
  
  if (vaf_noise) {
    
    for (i in seq_along(s$samples)) {
      
      # get arguments of sample
      smp = s$samples[[i]]
      
      stopifnot(all(c("min_vaf","min_reads","depth_model","depth") %in% names(smp)))
      min_vaf = smp$min_vaf
      min_reads =  smp$min_reads
      depth_model = smp$depth_model
      depth = smp$depth
      
      stopifnot(depth_model == 1)
      
      # depth to test
      dp_vals = seq(0, 1e4) # long enough?
      min_reads_dp = sapply(dp_vals, function(x) max(c(min_reads, ceiling(x*min_vaf))))
      
      p_dp = dpois(dp_vals, depth) # prop to see given depth
      stopifnot(all.equal(sum(p_dp), 1.0))
      p_detect = pbinom(min_reads_dp, dp_vals, vaf_noise, lower.tail=FALSE) # prop to see sufficent noise at given depth
      
      n_per_dp_exp = genome_size * p_dp * p_detect
      n_per_dp_smp = rpois(length(n_per_dp_exp), n_per_dp_exp)
      
      if (sum(n_per_dp_smp) > 0) {
        
        dp = unlist(lapply(seq_along(n_per_dp_smp), function(i) rep(dp_vals[i], n_per_dp_smp[i])))
        dp_min = unlist(lapply(seq_along(n_per_dp_smp), function(i) rep(min_reads_dp[i], n_per_dp_smp[i])))
        alt = extraDistr::rtbinom(length(dp), dp, vaf_noise, a = dp_min)
        
        data_noise =
          data.frame(
            clone = rep(NA, length(alt)),
            alt = alt,
            depth = dp,
            id = paste0("N", seq(n_noisy_sites, n_noisy_sites + length(alt)-1)),
            ccf = rep(0.0, length(alt)),
            vaf = alt / dp,
            cn = rep(2, length(alt))
          )
        
        n_noisy_sites = n_noisy_sites + length(alt)
        
        s$samples[[i]]$mutation_data =
          rbind(s$samples[[i]]$mutation_data, data_noise)
        
      }
    }
  }
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check detection status of all samples
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  for (i in seq_along(s$samples)) {
    
    # get samples argument to determine if variant was detected
    smp = sample_list[[i]]
    if (!"min_vaf" %in% names(smp)) { smp$min_vaf = 0.05 }
    if (!"min_reads" %in% names(smp)) { smp$min_reads = 2 }
    
    stopifnot(all(c("min_vaf","min_reads") %in% names(smp)))
    min_vaf = smp$min_vaf
    min_reads =  smp$min_reads
    
    # add detection status to variants
    s$samples[[i]]$mutation_data =
      s$samples[[i]]$mutation_data %>% 
      dplyr::mutate(detected = (vaf >= min_vaf | min_vaf == 0) & alt >= min_reads)
    
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Keep or add detected variants
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # keep only detected variants
  detected_variants = 
    s$samples %>% 
    lapply(function(x) x$mutation_data[x$mutation_data$detected,]) %>% 
    do.call(what=rbind) %>% 
    dplyr::select(clone, id, cn) %>% 
    unique()
  
  for (i in seq_along(s$samples)) {
    
    smp = s$samples[[i]]
    s_data =  smp$mutation_data
    
    stopifnot(all(c("depth_model","depth") %in% names(smp)))
    depth_model = smp$depth_model
    depth =  smp$depth
    
    if (depth_model == 1) {
      get_dp = function(n) rpois(n, depth) # poisson 
    } else if (depth_model == 1) {
      if (NROW(s_data) < 50)
        warning("Few variants to resample for depth in 'add_all_samples'.\n")
      if (NROW(s_data) < 5)
        stop("To few variants to resample for depth in 'add_all_samples'.\n")
      get_dp = function(n) sample(s_data$depth, n, replace=1) # resample data
    } else if (depth_model == 3) {
      get_dp = function(n) rep(depth, n) # constant
    } else {
      stop("Unknown depth model.\n")
    }
    
    # add missing variants
    missing_data_to_add = 
      detected_variants[ !detected_variants$id %in% s_data$id,] %>% 
      dplyr::mutate(alt=0, vaf=0.0, ccf=0, detected=TRUE) %>% 
      dplyr::mutate(depth=get_dp(length(alt))) %>% 
      dplyr::select(clone, alt, depth, id, ccf, vaf, cn)
    
    
    # drop undetected variants not found in another sample
    s$samples[[i]]$mutation_data =
      s$samples[[i]]$mutation_data %>% 
      dplyr::filter(id %in% detected_variants$id | detected) %>% 
      dplyr::select(-detected)
    
    
    s$samples[[i]]$mutation_data = 
      rbind(s$samples[[i]]$mutation_data, missing_data_to_add)
    
  }
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Optionally add noise to VAF data
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  if (vaf_noise) {
    
    for (i in seq_along(s$samples)) {
      
      add_noise = grepl("^N", s$samples[[i]]$mutation_data$id) & s$samples[[i]]$mutation_data$vaf > 0
      dp = s$samples[[i]]$mutation_data$depth
      n_alt = s$samples[[i]]$mutation_data$alt
      n_alt_back = rbinom(length(n_alt), n_alt, vaf_noise / 4)
      
      n_ref = dp - n_alt - n_alt_back
      n_to_alt = rbinom(length(n_ref), n_ref, vaf_noise)
      
      new_alt =  n_alt - n_alt_back + n_to_alt
      new_alt[!add_noise] = n_alt[!add_noise]
      s$samples[[i]]$mutation_data$alt = new_alt
      s$samples[[i]]$mutation_data$vaf = new_alt / dp
    }
    
  }
  
  return(s)
}



#' #' Title
#' #'
#' #' @param s 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' add_cn_states = function(s, , seed=123) {
#'   
#'   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#'   # Alternative config ####
#'   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#'   
#'   if (isTRUE(f == "WES")) {
#'     f = 0.011 # 1.1%
#'   }
#'   
#'   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#'   # Check arguments
#'   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#'   
#'   al = checkmate::makeAssertCollection()
#'   if (is.data.frame(s)) {
#'     checkmate::assertTRUE("id" %in% colnames(s), add=al)
#'   } else {
#'     checkmate::assertClass(s, "CHESS_S3", add=al)
#'   }
#'   checkmate::assertNumeric(f, lower=0, upper=1, len=1, any.missing=FALSE, add=al)
#'   checkmate::reportAssertions(al)
#'   
#'   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#'   # Subset mutation data ####
#'   # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#'   
#'   w = c(f, 1-f) # weight
#'   
#'   if (is.data.frame(s)) {
#'     
#'     # data.frame
#'     s = s %>% dplyr::filter(mutation_id_to_state(id, weight=w, seed=seed) == 1)
#'     
#'   } else {
#'     
#'     # CHESS S3 object
#'     for (i in seq_along(s$samples)) {
#'       s$samples[[i]]$mutation_data = 
#'         s$samples[[i]]$mutation_data %>% 
#'         dplyr::filter(mutation_id_to_state(id, weight=w, seed=seed) == 1)
#'       
#'       s$samples[[i]]$number_mutations = 
#'         NROW(s$samples[[i]]$mutation_data)
#'     }
#'     
#'   }
#'   
#'   
#' }


#' Random sampling of mutations
#'
#' @param s CHESS S3 object.
#' @param f Numeric value defining fraction of variants to keep (0<f<=1). Alternative options: "WES": whole-exome sequencing, i.e. 1.1%.
#' @param seed (optional) seed to use for random subsetting (default: 123).
#'
#' @return
#' @export
#'
#' @examples
subset_mutations = function(s, f, seed=123) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Alternative config ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  if (isTRUE(f == "WES")) {
    f = 0.011 # 1.1%
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  if (is.data.frame(s)) {
    checkmate::assertTRUE("id" %in% colnames(s), add=al)
  } else {
    checkmate::assertClass(s, "CHESS_S3", add=al)
  }
  checkmate::assertNumeric(f, lower=0, upper=1, len=1, any.missing=FALSE, add=al)
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Subset mutation data ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  w = c(f, 1-f) # weight
  
  if (is.data.frame(s)) {
    
    # data.frame
    s = s %>% dplyr::filter(mutation_id_to_state(id, weight=w, seed=seed) == 1)
    
  } else {
    
    # CHESS S3 object
    for (i in seq_along(s$samples)) {
      s$samples[[i]]$mutation_data = 
        s$samples[[i]]$mutation_data %>% 
        dplyr::filter(mutation_id_to_state(id, weight=w, seed=seed) == 1)
      
      s$samples[[i]]$number_mutations = 
        NROW(s$samples[[i]]$mutation_data)
    }
    
  }
  
  return(s)
}


#' Internal function to map mutations to a state vector
#'
#' @param id Mutation id
#' @param weight Vector of number weights
#' @param seed optional seed to use for random permutation of mapping.
#'
#' @return vector with length equal to id with state index
#'
#' @examples table(mutation_id_to_state(as.character(1:1000), c(4, 1, 5)))
mutation_id_to_state = function(id, weight, seed=123) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  checkmate::assertCharacter(id, add=al)
  checkmate::assertNumeric(weight, lower = 0, min.len = 1, finite = TRUE, any.missing = FALSE, add=al)
  checkmate::assertNumeric(seed, any.missing = FALSE, null.ok = FALSE, finite = TRUE, len = 1, add=al)
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  if (is.null(id)) return(NULL)
  
  mut_q = # mutation quantile (hash of variant id, scaled to (0,1])
    as.character(id) %>% 
    sapply(digest::digest, algo="xxhash32", raw=TRUE, seed=seed) %>% 
    (function(x) as.numeric(paste0("0x", x)) / (2^32))
  
  # normalize weights and get cumsum
  weight = weight / sum(weight) 
  weight_cs = cumsum(weight)
  
  # mutation to state
  state_of_mut = sapply(mut_q, function(x) which(x <= weight_cs)[1])
  names(state_of_mut) = as.character(id)
  
  return(state_of_mut) 
}

