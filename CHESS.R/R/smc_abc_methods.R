#' Calculate particle weights
#'
#' @param d A nmatrix (2D) or array (3D) of distance values.
#' @param epsilon The distance threshold of accepted particles.
#'
#' @return Particle weights.
#' @export
#'
calc_w =
  function(d, epsilon) {
    ac = checkmate::makeAssertCollection()
    checkmate::assertArray(d, mode = "numeric", min.d=2, max.d=3, add=ac)
    checkmate::assertNumeric(d, lower=0, add=ac)
    checkmate::assertNumber(epsilon, lower=0, add=ac)
    checkmate::reportAssertions(ac)
    
    below_eps = d < epsilon
    
    if (length(dim(d)) == 3) {
      below_eps = apply(below_eps, c(1,2), sum, na.rm=TRUE)
    }
    
    w = apply(below_eps, 1, sum, na.rm=TRUE)
    w / sum(w)
  }




#' Calculate  effective sample size (ESS)
#'
#' @param weight A numeric vector of particle weights. This will be normalized. 
#'
#' @return The effective sample size of the current particle population.
#' @export
#'
calc_ess =
  function(weight) {
    checkmate::assertNumeric(weight, lower=0, any.missing=FALSE, finite=TRUE)
    weight = weight / sum(weight)
    sum(weight ^ 2) ^ -1
  }


#' Calculate distance matrix of a particle set
#'
#' @param p The particle set.
#' @param y The target tree.
#' @param rho The distance function. Must return a proper distance. 
#' @param multi_summary (optional) Function to simplified multiple observations (default: none).
#' @param use_precalc_d (optional) Logical indicating if pre-calculated distance should be trusted (default: true).
#' @param simplify (optional) Logical indicating if the distance array should be simplified to a matrix (default: true).
#'
#' @return A 2d matrix (particles, realizations * multi_summary(observations)) or a 3d array (particles, realizations, observations) if "simplify" is FALSE and "multi_summary" is NULL. Missing values (i.e. NA) will be set to Inf!
#' 
d_matrix_from_particles_pb =
  function(p, y, rho, multi_summary=NULL, use_precalc_d=TRUE, simplify=TRUE, n_cores=1) {
    
    # check arguments
    ac = checkmate::makeAssertCollection()
    checkmate::assertClass(y, "phylo", null.ok = TRUE, add=ac)
    checkmate::assertFunction(rho, null.ok = TRUE, nargs=2, add=ac)
    checkmate::assertFunction(multi_summary, nargs=1, null.ok = TRUE, add=ac)
    checkmate::assertFlag(use_precalc_d, add=ac)
    checkmate::assertFlag(simplify,  add=ac)
    checkmate::reportAssertions(ac)
    
    # *quick* check of distance metric (not exhaustive, but worth to verify) 
    if (!is.null(rho)) {
      t1 = ape::rtree(5); t1$tip.label[1] = "GL"
      t2 = ape::rtree(5); t2$tip.label = t1$tip.label
      checkmate::assertTRUE(all.equal(rho(t1, t2), rho(t2, t1))) # test that symetric
      checkmate::assertTRUE(rho(t1, t1) == 0.0)         # test that indiscernible
      checkmate::assertTRUE(rho(t1, t2) > 0.0)          # test that positive
    }
    
    # *quick* check of multi_summary function
    if (!is.null(multi_summary)) {
      checkmate::assertNumber(multi_summary(c(1:10)))
    }
    
    # wrapper around distance function:
    rho_ = # modified rho catching missing trees
      function(x, y) { 
        if (is.null(x)) return(NA)
        if (length(x) == 1 & is.na(x[1])) return(NA)
        if (is.null(x$tree)) return(NA)
        rho(x$tree, y)
      }
    
    
    # Calculate the distance array
    dmat = 
      abind::abind(along=0, pbmcapply::pbmclapply(p, function(a) {
        if ("distances" %in% names(a) & use_precalc_d) { # use pre-calculated distances
          abind::abind(along=1, lapply(a$distances, function(b) {
            abind::abind(along=2, as.list(b))
          }))
        } else { # re-calculated distances
          abind::abind(along=1, lapply(a$observations, function(b) {
            if (is.list(b) & !"tree" %in% names(b)) {
              return(abind::abind(along=2, lapply(b, function(c) { rho_(c, y) })))
            } else {
              return(rho_(b, y))
            }
          }))
        }
      }, mc.cores=n_cores))

   # Summaries multiple observations of one realization:
    if (length(dim(dmat)) == 3 & simplify) { # 3d array
      if (!is.null(multi_summary)) { # use summary function
        dmat = apply(dmat, c(1,2), multi_summary)
      } else { # bind along second dimension
        dmat_bound = NULL
        for (i in seq_len(dim(dmat)[3])) {
          dmat_bound = cbind(dmat_bound, dmat[,,i])
        }
        dmat = dmat_bound
      }
    }
    dmat[is.na(dmat)] = Inf
    
    return(dmat)
  }


#' Calculate distance matrix of a particle set
#'
#' @param p The particle set.
#' @param y The target tree.
#' @param rho The distance function. Must return a proper distance. 
#' @param multi_summary (optional) Function to simplified multiple observations (default: none).
#' @param use_precalc_d (optional) Logical indicating if pre-calculated distance should be trusted (default: true).
#' @param simplify (optional) Logical indicating if the distance array should be simplified to a matrix (default: true).
#'
#' @return A 2d matrix (particles, realizations * multi_summary(observations)) or a 3d array (particles, realizations, observations) if "simplify" is FALSE and "multi_summary" is NULL. Missing values (i.e. NA) will be set to Inf!
#' 
d_matrix_from_particles =
  function(p, y, rho, multi_summary=NULL, use_precalc_d=TRUE, simplify=TRUE, n_cores=1) {
    
    # check arguments
    ac = checkmate::makeAssertCollection()
    checkmate::assertClass(y, "phylo", null.ok = TRUE, add=ac)
    checkmate::assertFunction(rho, null.ok = TRUE, nargs=2, add=ac)
    checkmate::assertFunction(multi_summary, nargs=1, null.ok = TRUE, add=ac)
    checkmate::assertFlag(use_precalc_d, add=ac)
    checkmate::assertFlag(simplify,  add=ac)
    checkmate::reportAssertions(ac)
    
    # *quick* check of distance metric (not exhaustive, but worth to verify) 
    if (!is.null(rho)) {
      t1 = ape::rtree(5); t1$tip.label[1] = "GL"
      t2 = ape::rtree(5); t2$tip.label = t1$tip.label
      checkmate::assertTRUE(all.equal(rho(t1, t2), rho(t2, t1))) # test that symetric
      checkmate::assertTRUE(rho(t1, t1) == 0.0)         # test that indiscernible
      checkmate::assertTRUE(rho(t1, t2) > 0.0)          # test that positive
    }
    
    # *quick* check of multi_summary function
    if (!is.null(multi_summary)) {
      checkmate::assertNumber(multi_summary(c(1:10)))
    }
    
    # wrapper around distance function:
    rho_ = # modified rho catching missing trees
      function(x, y) { 
        if (is.null(x)) return(NA)
        if (length(x) == 1 & is.na(x[1])) return(NA)
        if (is.null(x$tree)) return(NA)
        rho(x$tree, y)
      }
    
    
    # Calculate the distance array
    dmat = 
      abind::abind(along=0, parallel::mclapply(p, function(a) {
        if ("distances" %in% names(a) & use_precalc_d) { # use pre-calculated distances
          abind::abind(along=1, lapply(a$distances, function(b) {
            abind::abind(along=2, as.list(b))
          }))
        } else { # re-calculated distances
          abind::abind(along=1, lapply(a$observations, function(b) {
            if (is.list(b) & !"tree" %in% names(b)) {
              return(abind::abind(along=2, lapply(b, function(c) { rho_(c, y) })))
            } else {
              return(rho_(b, y))
            }
          }))
        }
      }, mc.cores=n_cores))

   # Summaries multiple observations of one realization:
    if (length(dim(dmat)) == 3 & simplify) { # 3d array
      if (!is.null(multi_summary)) { # use summary function
        dmat = apply(dmat, c(1,2), multi_summary)
      } else { # bind along second dimension
        dmat_bound = NULL
        for (i in seq_len(dim(dmat)[3])) {
          dmat_bound = cbind(dmat_bound, dmat[,,i])
        }
        dmat = dmat_bound
      }
    }
    dmat[is.na(dmat)] = Inf
    
    return(dmat)
  }



#' SMC-ABC on a phylogenetic tree
#'
#' @param const_params A set of constant parameters for the function call of 'new_simulation' (see its manual for details: new_simulation?)
#' @param p_gen Particle generator function. These provide additional parameters for the function call of 'new_simulation' (this would be the 'param.matrix' and 'seed'). 
#' @param s_gen Sample generator function. This returns a sampling schema for a given simulation (first argument). 
#' @param target_y The target tree (object of class phylo.).
#' @param rho 
#' @param N 
#' @param M 
#' @param N_T 
#' @param min_acceptance_rate 
#' @param alpha 
#' @param marginal_kernel 
#' @param epsilon_min 
#' @param max_steps 
#' @param n_cores 
#' @param upper_bound 
#' @param lower_bound 
#' @param output_dir 
#' @param sample_tree 
#' @param epsilon 
#' @param M_sim 
#' @param min_distance_per_sim 
#' @param reobserve_resampled 
#' @param particle_collection_dir 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
smc_abc_tree =  function(const_params, p_gen, s_gen, target_y, rho, N=100, M=5, N_T=round(N/2), min_acceptance_rate=0, alpha=0.9, marginal_kernel=FALSE, epsilon_min=0, max_steps = 20, n_cores = (detectCores() - 2), upper_bound=NULL, lower_bound=NULL, output_dir=NULL, sample_tree=FALSE, M_sim=1, min_distance_per_sim=TRUE, reobserve_resampled=FALSE, particle_collection_dir=NULL, keep_closest_tree_only=TRUE, allow_increase_in_eps=FALSE, delta_epsilon_wait=5, delta_epsilon_min=0, simulation_manipulator=NULL, tree_manipulator=NULL, previous_state=NULL, n_cores_dist=n_cores, ...) {
  
  divider_text="\n=========================================================\n\n"
  
  # 
  # opt_file = file.path(output_dir, "options.rds")
  # if (file.exists(opt_file)) {
  # 
  #   opts = readRDS(opt_file)
  #   
  #   for (var in names(opts)) {
  #     assign(var, opts[[var]])
  #   }
  # }
  
  rho = tryCatch(get(rho), error = function(e) return(rho)) # try getting object with name in rho (i.e. get the named function)
  
  # check input
  stopifnot(all(is.numeric(c(N,M_sim,M,alpha,epsilon_min))))
  stopifnot(all(sapply(list(N,M_sim,M,alpha,epsilon_min,rho), length) == 1))
  stopifnot(all(sapply(list(rho, p_gen, s_gen), is.function)))
  stopifnot(!file.exists(output_dir) | length(list.files(output_dir)) == 0)
  dir.create(output_dir, FALSE, TRUE) # create dir.
  stopifnot(all(sapply(list(rho, p_gen, s_gen), is.function)))
  stopifnot(class(target_y) == "phylo")
  
  # check that s_gen accepts simulation as argument.
  # check that p_gen returns particles variable parameters.
  # check content of const_params returns particles variable parameters.
  # check content of p_gen.
  # check content of s_gen (e.g. does it return all samples in target_y?). 
  # check that rho works. 
  
  # test rho
  rho_test_res = 
    tryCatch({ rho(target_y, target_y) },
             error = function(e) {
               stop("Function 'rho' does not work correctly. See ?smc_abc_tree for details.")
             })
  
  if (!is.numeric(rho_test_res)) stop("Function 'rho' does not return a numeric distance.")
  if (length(rho_test_res) != 1) stop("Function 'rho' does not return a numeric distance.")
  if (rho_test_res != 0) stop("Function 'rho' does not return a numeric distance.")
  
  
  # save all local variables to element
  vars = ls()
  var_list = list()
  for (var in vars) {
    var_list[[var]] = get(var)
  }
  saveRDS(var_list, file.path(output_dir, paste0("options.rds")), version = 2)
  
  # header message
  cat("+=======================================================+\n")
  cat("|                        SMC ABC                        |\n")
  cat("+=======================================================+\n")
  cat("\n")
  cat(" Options:\n")
  cat("  - N =", N, "particles.\n")
  cat("  - M =", M, "realisation(s) per particle.\n")
  cat("  - M_sim =", M_sim, " observation(s) per realisation.\n")
  cat("  - N_T =", N_T, "minimum effective sample size.\n")
  cat("  - alpha =", alpha, "ESS reduction per step.\n")
  cat("  - delta_epsilon_min =", delta_epsilon_min, "target distance for", delta_epsilon_wait, "steps (termination condition).\n")
  cat("  - epsilon_min =", epsilon_min, "target distance (termination condition).\n")
  cat("  - max_steps =", max_steps, " maximum steps  (termination condition). \n")
  cat("  - min_acceptance_rate =", min_acceptance_rate, " (termination condition).\n")
  cat("  - n_cores =", n_cores, "\n")
  cat("  - output_dir =", output_dir, "\n")
  cat(divider_text)
  start_time = Sys.time()
  
  
  #####################
  # SMC variables
  #####################
  
  # params for smc
  if (!is.null(previous_state)) {
    
    wh_params_needed = c("n_step","avg_acceptance","delta_d","epsilon","rng_state","particles")
    wh_missing = !wh_params_needed %in% names(previous_state)
    
    if (any(wh_missing)) {
      stop(paste0("Missing one or several elements in the argument 'previous_state':\n", 
                  "  => ", paste0(wh_params_needed[wh_missing], collapse=", ")))
    }
      
    for (param in wh_params_needed) {
      assign(param, previous_state[[param]])
    }
    
    if (!is.null(rng_state)) {
      .Random.seed  = rng_state
    }

  } else {
    
    # collect initial particles
    particles = list()
    
    while (length(particles) != N) {
      
      particles_new =
        sample_particles( # draw particles from prior
          constant_params = const_params, 
          generator_sim_params = p_gen,
          N = N - length(particles)
        )
      
      if (!is.null(lower_bound) | !is.null(upper_bound)) {
        
        keep = rep(TRUE, length(particles_new))
        
        for (i in seq_along(particles_new)) {
          
          params_i =  particles_new[[i]]$param.matrix
          
          # check particle above lower bound
          if (!is.null(lower_bound)) {
            
            stopifnot(all.equal(dim(params_i), dim(lower_bound)))
            stopifnot(all.equal(dimnames(params_i), dimnames(lower_bound)))
            
            if (any(params_i < lower_bound, na.rm = TRUE)) {
              keep[i] = FALSE
            }
          }
          
          # check particle below upper bounds
          if (!is.null(upper_bound)) {
            
            stopifnot(all.equal(dim(params_i), dim(upper_bound)))
            stopifnot(all.equal(dimnames(params_i), dimnames(upper_bound)))
            
            if (all(params_i > upper_bound, na.rm = TRUE)) {
              keep[i] = FALSE
            }
          }
          
        }
        
        if (sum(!keep) == N) stop("All particles reject. Bad prior?")
        particles_new = particles_new[keep]
      } 
      
      particles = c(particles, particles_new)
      
    }
    
    # add M observations to particles
    if (!"observations" %in% names(particles[[1]])) {
      
      cat("Creating initial particle population:\n\n")
      
      particles = 
        observe_particles( # create M x M_sim realizations of particles
          particles = particles,
          generator_sample_params = s_gen, 
          M = M, 
          n_cores = n_cores, 
          sample_tree = sample_tree,
          M_sim = M_sim, 
          rho=rho, 
          target_tree = target_y,
          dump_dir = particle_collection_dir, 
          keep_closest_tree_only = keep_closest_tree_only,
          simulation_manipulator = simulation_manipulator,
          tree_manipulator = tree_manipulator,
          ...
        )
    }
    
    n_step = 0
    avg_acceptance = 1
    delta_d = c(Inf)
    epsilon = Inf
    

  }

  d_matrix = d_matrix_from_particles(particles, target_y, rho, n_cores = n_cores_dist)
  W = calc_w(d_matrix, epsilon)
  ESS = calc_ess(W)
  
  
  # parameter labels
  id_row = rownames(particles[[1]]$param.matrix)
  id_col = seq_len(NCOL(particles[[1]]$param.matrix))
  id_param = paste0(rep(id_row, length(id_col)), rep(id_col,each=length(id_row)))
  
  
  # find variable parameters
  paramter_array =  particles %>% lapply("[[", "param.matrix") %>% abind::abind(along = 3) 
  
  idx_var_params =  # index of varied parameters, some are constants.
    paramter_array %>%
    apply(c(1,2), function(x) length(unique(x)) != 1) %>% 
    c() %>% 
    magrittr::set_names(id_param)
  
  tied_params = # e.g. always equal
    paramter_array %>%
    apply(c(1,3), function(x) length(unique(x)) == 1) %>% 
    apply(1, all) 
  
  # index of tied parameters:
  tied_params_idx = seq_along(tied_params) #%>% magrittr::set_names(names(tied_params))
  tied_params_idx[!tied_params] = NA
  tied_params_idx = rep(tied_params_idx, dim(paramter_array)[2])
  names(tied_params_idx) = id_param
  
  # unvaried parameters are not tied:
  tied_params_idx[!idx_var_params] = NA
  
  #  tied parameters are not varied:
  idx_var_params[which(seq_along(tied_params_idx) != tied_params_idx)] = FALSE
  
  # ...
  final_parameter_index = idx_var_params
  final_parameter_index[final_parameter_index] = seq_len(sum(final_parameter_index))
  wh = !is.na(tied_params_idx)
  final_parameter_index[wh] = final_parameter_index[tied_params_idx[wh]]
  
  #test = names(which(idx_var_params))[final_parameter_index] #%>% set_names(names())
  #names(test) = names(which(final_parameter_index != 0))
  
  # lower bound
  if (is.null(lower_bound)) { # default lower bound is 0
    lower_bound = rep(0, sum(idx_var_params))
  } else {
    lower_bound = lower_bound[idx_var_params]
  }
  names(lower_bound) = id_param[idx_var_params]
  
  # upper bound
  if (is.null(upper_bound)) { # default upper bound is Inf
    upper_bound = rep(Inf, sum(idx_var_params))
  } else {
    upper_bound = upper_bound[idx_var_params]
  }
  names(upper_bound) = id_param[idx_var_params]
  
  
  #####################
  # SMC implementation
  #####################
  
  dir.create(file.path(output_dir, "rejected"), FALSE, TRUE)
  
  if (is.null(previous_state)) {
    
    # store chain of particles and parameters
    state =
      list(
        n_step = n_step,
        particles = particles,
        epsilon = epsilon,
        weight = W,
        ess = ESS,
        acceptance_rate = avg_acceptance,
        rng_state = .Random.seed,
        delta_d = delta_d
      )
    
    
    # optionally, save output to file
    if (!is.null(output_dir)) {
      out_file = file.path(output_dir, paste0("state_step_", n_step,".rds"))
      saveRDS(state, out_file, version = 2)
    } else {
      chain = list()
      rejected_particles = list()
      chain = c(chain, list(state))
    }
  }
  
  
  
  while (TRUE) { # do smc steps until epsilon reached
    
    # termination criteria
    n_step = n_step + 1 
    if (is.infinite(N_T)) break
    if (n_step >= max_steps) N_T = Inf # trigger one more resample step.
    if (avg_acceptance <= min_acceptance_rate) N_T = Inf
    if (epsilon <= epsilon_min) N_T = Inf
    if (max(tail(delta_d, delta_epsilon_wait), na.rm=TRUE) < delta_epsilon_min) N_T = Inf
    if (is.infinite(N_T)) break

    cat("Step ", n_step, ":\n\n", sep="")
    cat(" *** Updating epsilon ***\n")
    cat("  => last rel. diff ESS: ")
    cat(paste0(signif(tail(delta_d, delta_epsilon_wait)*100, 3), "%", collapse=", "))
    cat("\n")
    
    # find new distance limit, epsilon
    if (ESS >= N_T) { # if the ESS is above the minimum sampling size
      d_max = max(d_matrix[!is.infinite(d_matrix)])
      fun = function(e) {val = calc_ess(calc_w(d_matrix, e)) - ESS * alpha; if (val >= 0) val + ESS else val}
      #epsilon_n = optimize(fun, lower=min(d_matrix), upper=d_max)$minimum
      epsilon_n = uniroot(fun, c(sort(unique(d_matrix))[2], upper=d_max))$root
      epsilon_n = min(c(epsilon_n, epsilon), na.rm=TRUE)
    } else { # do not reduce distance even further if the ESS is already below the minimum
      epsilon_n = epsilon
      #cat("  => ESS < N_T (", ESS, " < ", N_T, ") : not reducing epsilon any further.\n", sep="")
      cat("  => ESS < N_T: not reducing epsilon any further.\n", sep="")
    }
    new_ess = calc_ess(calc_w(d_matrix, epsilon_n))
    
    
    
    # updated parameters
    d_eps = epsilon - epsilon_n
    delta_d = c(delta_d, d_eps / epsilon)
    d_eps_rel = signif(d_eps / epsilon * 100, 2)
    epsilon = epsilon_n
    W = calc_w(d_matrix, epsilon)
    ess_n = calc_ess(W)
    d_ess_rel = signif(100 - ess_n / ESS * 100, 2) 
    ESS = ess_n
    
    n_failed = sum(is.na(d_matrix))
    if (n_failed) {
      frac_failed = signif(n_failed / prod(dim(d_matrix)) * 100, 2)
      cat("  => Some particles failed observation: N = ", n_failed, " (", frac_failed, "%)\n",sep="")
    }
    cat("  => Epsilon: ",epsilon," (delta: ",d_eps,", ",d_eps_rel,"%)\n",sep="")
    cat("  => ESS: ", ESS, " (delta: ", d_ess_rel, "%)\n", sep="")
    if (ESS < N_T) 
      cat("  => ESS < N_T (", round(ESS), " < ", N_T, ")\n", sep="")
    cat("\n")
    
    
    #####################
    # resampling step?
    #####################
    
    if (ESS < N_T) { # resampling step
      cat(" *** Resampling ***\n")
      
      ESS_pre_sampling = ESS
      idx =  sample(seq_along(W), replace=TRUE, prob=W)
      
      if (reobserve_resampled) {
        
        idx_dup = duplicated(idx)
        particles = particles[idx]
        
        particles[idx_dup] = 
          observe_particles( # create M x M_sim realizations of particles
            particles = particles[idx_dup],
            generator_sample_params = s_gen, 
            M = M, 
            n_cores = n_cores, 
            sample_tree = sample_tree,
            M_sim = M_sim, 
            rho=rho, 
            target_tree = target_y,
            dump_dir = particle_collection_dir, 
            keep_closest_tree_only=keep_closest_tree_only,
            simulation_manipulator=simulation_manipulator,
            tree_manipulator = tree_manipulator,
            ...
          )
        
        d_matrix[idx_dup,] = d_matrix_from_particles(particles[idx_dup], target_y, rho, n_cores = n_cores_dist)
        
      } else {
        particles = particles[idx]
        d_matrix = d_matrix[idx,]
      }
      
      
      # report new ess after pretubation
      ess_pp = calc_ess(calc_w(d_matrix, epsilon))
      d_ess_rel = signif((ess_pp - ESS) / ESS * 100, 2)
      cat("  => ESS: ", ess_pp, " (delta: ", d_ess_rel, "%)\n", sep="")
      
      # update weight and ESS
      W = rep(1 / length(W), length(W))
      ESS = N
    }
    
    #####################
    # pertubation step
    #####################
    
    cat(" *** Pertubation ***\n")
    
    pertube = W != 0 # don't pertube particles with W == 0.
    cat("  => Proposing new parameters for", sum(pertube), " particles:\n\n")
    
    theta = 
      particles %>% 
      lapply(function(x) c(x$param.matrix)[idx_var_params]) %>% 
      do.call(what=rbind) %>% 
      magrittr::set_colnames(names(idx_var_params[idx_var_params]))
    
    # mutation rate from tree scaling
    est_mutation_rate = unlist(sapply(particles, function(x) ifelse(is.null(x$mutation_rate), NA, x$mutation_rate)))
    if (all(is.na(est_mutation_rate))) {
      theta_mod = theta
    } else {
      theta_mod = cbind(theta, mutation_rate=est_mutation_rate)
    }
    
    # estimate covariance
    wh_use = apply(!is.na(theta_mod), 1, all)
    cov_est = cov.wt(theta_mod[wh_use,,drop=FALSE], W[wh_use], cor = TRUE)
    cov = cov.wt(theta[wh_use,,drop=FALSE], W[wh_use], cor = TRUE)
    print_digits = 2
    
    cat("  - median_theta:")
    print(knitr::kable(signif(apply(theta_mod, 2, median, na.rm=TRUE), print_digits)))
    cat("\n\n")
    
    cat("  - mu_theta:")
    print(knitr::kable(signif(apply(theta_mod, 2, weighted.mean, w=W), print_digits)))
    cat("\n\n")
    
    cat("  - cov_theta:")
    print(knitr::kable(signif(cov_est$cov, print_digits)))
    cat("\n\n")
    
    cat("  - cor_theta:")
    print(knitr::kable(signif(cov_est$cor, print_digits)))
    cat("\n\n")
    
    
    # sigma for mvnorm
    sigma = (cov$cov + diag(ncol(cov$cov)) * 0.01) * 2
    if (marginal_kernel) {
        sigma[upper.tri(sigma) | lower.tri(sigma)] = 0
    }
    
    # propose new theta
    wh_pertube = which(pertube)
    theta_star = do.call(cbind, lapply(wh_pertube, function(i) {
      c(tmvtnorm::rtmvnorm(theta[i,],  n=1, sigma=sigma, 
                           lower=lower_bound, upper=upper_bound))
    }))
    
    q_ratios =  sapply(seq_along(wh_pertube), function(i) {
      q1 = tmvtnorm::dtmvnorm(theta[wh_pertube[i],], theta_star[,i], 
                              sigma=sigma, lower=lower_bound, upper=upper_bound)
      q2 = tmvtnorm::dtmvnorm(theta_star[,i], theta[wh_pertube[i],], 
                              sigma=sigma, lower=lower_bound, upper=upper_bound)
      q1 / q2
    })
    
   
    
    # update copy of particles with theta star
    particles_star = particles[pertube]
    for (i in seq_along(particles_star)) {
      wh = final_parameter_index != 0
      idx = final_parameter_index[wh]
      particles_star[[i]]$param.matrix[wh] = theta_star[idx,i]
    }
    
    # reorder subclones based on their respective start times
    for (i in seq_along(particles_star)) {
      params = particles_star[[i]]$param.matrix
      ord = order(params["clone_start_times",])
      particles_star[[i]]$param.matrix[,ord] = params[,ord]
    }
    
    # draw observations and get accepted set
    particles_star =
      observe_particles(
        particles = particles_star, 
        generator_sample_params = s_gen,
        M = M, 
        n_cores = n_cores, 
        sample_tree = sample_tree, 
        M_sim =  M_sim, 
        rho = rho, 
        target_tree = target_y, 
        dump_dir = particle_collection_dir, 
        keep_closest_tree_only = keep_closest_tree_only,
        simulation_manipulator = simulation_manipulator,
        tree_manipulator = tree_manipulator,
        ...
      )
    
    
    d_matrix_star = d_matrix_from_particles(particles_star, target_y, rho, n_cores = n_cores_dist)
    
    # mh ratio:
    n_star = apply(d_matrix_star < epsilon, 1, sum)
    n = apply(d_matrix[pertube,,drop=FALSE] < epsilon, 1, sum)
    p_ratio =  n_star / n
    p_accept = p_ratio * q_ratios; 
    p_accept[p_accept > 1] = 1; p_accept[is.nan(p_accept)] = 0
    #pi_star / pi? 
    
    # reject or accept proposals
    accept_star = as.logical(rbinom(length(n_star), 1, p_accept))
    particles[pertube][accept_star] = particles_star[accept_star]
    d_matrix[which(pertube)[accept_star],] = d_matrix_star[accept_star,]
    pertube[pertube] = !accept_star
    
    # store rejected particles
    if (!is.null(output_dir)) { # optionally, save output to file
      out_file = file.path(output_dir, "rejected", paste0("S", n_step,".rds"))
      saveRDS(particles_star[!accept_star], out_file, version = 2)
    } else {
      rejected_particles = c(rejected_particles, list(particles_star[!accept_star]))
    }
    
    cat("   => Accepting ", sum(accept_star), "/",length(accept_star), 
        " (", round(mean(accept_star)*100), "%) proposal(s).\n", sep="")
    
    avg_acceptance = mean(accept_star)
    # update ess and weight after pertubation
    W = calc_w(d_matrix, epsilon)
    ess_pp = calc_ess(W)
    d_ess_rel = signif((ess_pp - ESS) / ESS * 100, 2)
    cat("  => ESS: ", ess_pp, " (delta: ", d_ess_rel, "%)\n\n", sep="")
    ESS = ess_pp
    cat(divider_text)
    
    #####################
    # save state
    #####################
    
    # store chain of particles and parameters
    state =
      list(
        n_step = n_step,
        particles = particles,
        epsilon = epsilon,
        weight = W,
        ess = ESS,
        acceptance_rate = avg_acceptance,
        rng_state = .Random.seed,
        delta_d = delta_d
      )
    
    # optionally, save output to file
    if (!is.null(output_dir)) {
      out_file = file.path(output_dir, paste0("state_step_", n_step,".rds"))
      saveRDS(state, out_file, version = 2)
    } else {
      chain = c(chain, list(state))
    }
    
  }
  
  # 
  rt_s = as.numeric(Sys.time()) - as.numeric(start_time)
  rt_h = round(rt_s / 3600, 1)
  cat("Finished. Total time spend: ", round(rt_s),"s (", rt_h, "h).\n", sep="")
  
  if (!is.null(output_dir)) {
    invisible(TRUE)
  } else {
    return(list(chain=chain, rejected_particles=rejected_particles))
  }
}




get_equivalent_single_cell_tree = function(s, setup, depth_cutoff=10, include_lowdepth=TRUE, include_node_labels=FALSE, insert_lp_above=FALSE) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  
  checkmate::assertClass(s, "Rcpp_CHESS_Universe_object", add=al)
  checkmate::assertList(setup, add=al) # more tests needed ...
  checkmate::assertNumeric(depth_cutoff, any.missing = FALSE, null.ok = FALSE, finite = TRUE, len=1, add=al)
  checkmate::assertFlag(include_lowdepth, add=al)
  checkmate::assertFlag(include_node_labels, add=al)
  
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # optionally remove low depth samples
  is_low_depth = sapply(setup, function(x) isTRUE(x$depth < depth_cutoff))
  names(is_low_depth) = names(setup)
  if (!include_lowdepth) {
    setup = setup[!is_low_depth]
    is_low_depth = is_low_depth[!is_low_depth]
  }
  
  # compare to single cell tree from phylogeny:
  x_cells = c()
  y_cells = c()
  z_cells = c()
  
  for (i in seq_along(setup)) {
    
    smp = setup[[i]]
    
    if (!"center" %in% names(smp) & "sample_coord" %in% names(smp)) { # alternative names
      tmp_coord = cbind(smp$sample_coord[1:3], smp$sample_coord[4:6])
      smp$center = floor(apply(tmp_coord, 1, mean))
      smp$diameter = apply(tmp_coord, 1, function(x) diff(range(x))) + 1
    }
    
    stopifnot(all(c("diameter","center") %in% names(smp)))
    stopifnot(length(smp$center) == 3 & length(smp$diameter) == 3)
    if (!all(smp$diameter == 1)) stop("Can't sample samples with d>1 ")
    
    smp$center = smp$center - 1 # 0 index
    if (s$CellType(smp$center[1], smp$center[2], smp$center[3]) == 0) next()
    x_cells = c(x_cells, smp$center[1])
    y_cells = c(y_cells, smp$center[2])
    z_cells = c(z_cells, smp$center[3])
    
  }
  coord_label = paste0("x", x_cells, "_y", y_cells, "_z", z_cells)
  stopifnot(!any(duplicated(coord_label)))
  
  # return null if tree construction is not possible
  if (length(x_cells) < 2) {
    return(NULL)
  }
  
  # get the tree
  if (include_node_labels) {
    single_cell_tree_data = s$SingleCellTreeWithNodeIds(x_cells, y_cells, z_cells)
    single_cell_tree_data_ = s$SingleCellTree(x_cells, y_cells, z_cells)
  } else {
    single_cell_tree_data = s$SingleCellTree(x_cells, y_cells, z_cells)
  }

  tree_object = ape::read.tree(text=single_cell_tree_data$tree_string)
  
  if (!include_node_labels) {
    # replace labels with the ones in the setup object
    tree_object = ape::collapse.singles(tree_object) # collapse all single nodes
    coord_label_tree = gsub(".*_type[[:digit:]]+_", "", tree_object$tip.label)
    tree_object$tip.label = names(setup)[match(coord_label_tree, coord_label)]
    tree_object$node.label = NULL
  } else {
    tree_object_ = ape::read.tree(text=single_cell_tree_data_$tree_string)
    stopifnot(tree_object_$node.label == tree_object$node.label)
    coord_label_tree = gsub(".*_type[[:digit:]]+_", "", tree_object_$tip.label)
    mapper_bc_to_id = tree_object$tip.label
    names(mapper_bc_to_id) = names(setup)[match(coord_label_tree, coord_label)]
    tree_object_ = ape::root.phylo(add.tips(tree_object_, "GL", getRoot(tree_object_), 0), "GL")
  }
  
  tree_object = ape::root.phylo(add.tips(tree_object, "GL", getRoot(tree_object), 0), "GL")
  
  # reassign lowpass samples
  if (any(is_low_depth)) {
    
    min_edge_length = floor(max(tree_object$edge.length)*0)
    label_dropped = names(which(is_low_depth))
    
    if (include_node_labels) {
      label_dropped = mapper_bc_to_id[label_dropped]
    }
    
    idx_dropped = which(tree_object$tip.label %in% label_dropped)
    
    for (i in idx_dropped) {
      
      c_edge = tree_object$edge[,2] == i # terminal edge 
      c_node = tree_object$edge[c_edge,1]
      
      wh_root = phangorn::getRoot(tree_object)
      label_dropped_this = label_dropped
      if (insert_lp_above) label_dropped_this = label_dropped[!label_dropped %in% tree_object$tip.label[i]]
      
      while (all(get_tips_below(c_node, tree_object) %in% label_dropped_this)) {
        c_edge = tree_object$edge[,2] == c_node # terminal edge 
        c_node = tree_object$edge[c_edge,1]
        if (sum(tree_object$edge[,2] == c_node)) break # fixes issues with root node?
      }
      
      new_label =  paste0(tree_object$tip.label[i], "added")
      stopifnot(!new_label %in% tree_object$tip.label)
      tree_object = phangorn::add.tips(tree_object, new_label, c_node, min_edge_length)
    }
    
    tree_object = ape::drop.tip(tree_object, label_dropped, collapse.singles = FALSE)
    tree_object$tip.label = gsub("added", "", tree_object$tip.label)
  }
  
  if (!include_node_labels){
    tree_object = ape::collapse.singles(tree_object)
  }
  
  return(tree_object) 
}


observe_particles =
  function(particles, generator_sample_params,  M, n_cores=(detectCores()-2), sample_tree=FALSE, M_sim=1, change_seed=TRUE, min_frac_clones=0.01, retries_simulation=3, retries_tree=5, rho=NULL, target_tree=NULL, dump_dir=NULL, keep_closest_tree_only=FALSE, simulation_manipulator=NULL, tree_manipulator = NULL, keep_trees=FALSE, ...) {
    
    if (M > 1 & !change_seed) stop("Multiple observations require a changed seed.\n")
    if (!change_seed) retries_simulation = 1
    if (M_sim>1 & is.null(generator_sample_params)) M_sim = 1
    stopifnot((is.null(rho) & is.null(target_tree)) | (!is.null(rho) & !is.null(target_tree)))
    
    
    # messages
    N = length(particles)
    cat("  => Collecting ", N, "x", M, "x", M_sim, " observations.\n", sep="")
    old_width = getOption("width")
    options(width=50)
    
    # indicator variables:
    precalculate_distances = !is.null(rho) & !is.null(target_tree)
    keep_min_obs_sim = precalculate_distances & M_sim > 1 & M > 1 & keep_closest_tree_only
    
    # helper functions
    if (sample_tree) {
      get_tree = 
        function(sim, loc, ...) {
          get_equivalent_single_cell_tree(sim$sim, loc)
        }
    } else {
      get_tree = 
        function(sim, loc, ...) {
          # add samples
          sim$samples = NULL
          sim = add_all_samples(sim, loc, ...)
          tryCatch(phylo_from_samples(sim, ...), error = function(e) return(NULL))
        }
    }
    
    
    if (precalculate_distances) {
      rho_ = function(x, y) { # modified rho catching missing trees
        tryCatch({
          if (is.null(x)) return(NA)
          if (length(x) == 1 & is.na(x[1])) return(NA)
          if (is.null(x$tree)) return(NA)
          rho(x$tree, y)
        }, error=function(e) return(NA))
      }
    }
    
    # 
    if (!is.null(target_tree)) {
      target_height = max(cophenetic(target_tree)["GL",])
    }
    
    # clear all particle observations
    for (i in seq_along(particles)) {
      particles[[i]]$observations = NULL
      particles[[i]]$distances = NULL
    }
    
    result = pbmcapply::pbmclapply(particles, function(particle) {
      
      particle$observations = list()
      if (precalculate_distances) particle$distances = list()
      
      for (k in seq_len(M)) { # simulate realisations of the same particle M times
        
        simulation = NULL
        
        tryCatch({  # try getting a acceptable simulation 'retries_simulation' times
          
          for (i in seq_len(retries_simulation)) { 
            
            if (change_seed) {
              particle$seed = round(runif(1, 0, 1e9))
            }
            
            # make sure some specific values are set
            particle_tmp = particle
            particle_tmp[c("record_history","verbose")] = list(FALSE, FALSE)
            
   
            # generate simulation, get clone stats
            simulation_try = do.call(CHESS::new_simulation, particle_tmp)
            n_clones = simulation_try$sim$CellCountsPerType
            frac_clones = n_clones / sum(n_clones)
            
            if (!is.null(dump_dir)) {
              tryCatch({ # store all simulations to a file
                stored_particle = particle
                stored_particle$frac_clones = n_clones
                wh_hash = !names(stored_particle) %in% c("observations","verbose")
                hash = digest::digest(stored_particle[wh_hash])
                timestamp = as.hexmode(round(as.numeric(Sys.time())))
                out_file = file.path(dump_dir, paste0(hash, ".", timestamp, "_sim.rds"))
                saveRDS(stored_particle, out_file, version = 2)
              }, error=function(e) return(NULL))
            }
            
            
            if (all(frac_clones > min_frac_clones) | (all(frac_clones[-1] > min_frac_clones) & length(frac_clones) > 2)) {
              simulation = simulation_try
              break
            }
          }
          
        }, error=function(e) {
          cat("Failed generation of simulation.\n")
          cat("Problematic particle:\n")
          particle_tmp$observations = NULL
          particle_tmp$distances = NULL
          dput(particle_tmp)
        })
        
        
        # optionally manipulate the simulation instance:
        if (!is.null(simulation_manipulator)) {
          simulation = simulation_manipulator(simulation)
        }

        
        if (is.null(simulation)) { 
          # failed to generate simulations,
          #  return nothing
          
          na_obs = list(seed=NULL, samples=NULL, tree=NULL, cell_counts=NULL)
          obs = rep(list(na_obs), M_sim)
          
        } else if (is.null(generator_sample_params)) { 
          # unobserved simulation,
          # return cell counts and seed only
          
          cell_count = simulation$sim$CellCountsPerType
          na_obs = list(seed=particle$seed, samples=NULL, tree=NULL, cell_counts=cell_count)
          obs = rep(list(na_obs), M_sim)
          
        } else { 
          # sucessful simulation,
          # observe trees and return full result list
          
          cell_count = simulation$sim$CellCountsPerType
          seed = particle$seed
          obs = rep(list(NULL), M_sim)
          
          get_obs = # helper function
            function(setup) {
              tryCatch({
                if (is.null(setup)) return(NULL)
                if ("setup" %in% names(setup)) smps = setup$setup else smps = setup
                tree = get_tree(simulation, smps, ...)
                if ("tree" %in% names(tree)) tree =tree$tree
                if (!all(names(setup$setup) %in% tree$tip.label)) { return(NULL) }
                if ("seed" %in% names(setup)) { setup$vars = NULL; setup$setup = NULL } # drop setup if they can be regenerated with the seed
                return(list(seed=seed, samples=setup, tree=tree, cell_counts=cell_count))
              }, error=function(e) return(NULL))
            }
          
          
          # try generating a valid tree up to 'retries_tree' times
          for (rep in seq_len(retries_tree)) { 
            wh_null = sapply(obs, is.null)
            n_null = sum(wh_null)
            if (n_null == 0) break
            sample_setups = replicate(n_null, tryCatch(generator_sample_params(simulation), error=function(e) return(NULL)), simplify=FALSE)
            obs[wh_null] = lapply(sample_setups, get_obs)
          } 
        
          # optionally, modify trees
          if (!is.null(tree_manipulator)) {
            for (i in seq_along(obs)) {
              if (!is.null(obs[[i]]$tree)) {
                obs[[i]]$tree = tree_manipulator(obs[[i]]$tree)
              }
            }
          }
          
          # failed observations
          na_obs = list(seed=seed, samples=NULL, tree=NULL, cell_counts=cell_count)
          wh_null = sapply(obs, is.null)
          obs[wh_null] = rep(list(na_obs), sum(wh_null))
         
          # scale generation tree to match target tree
          if (simulation$params$return_generation_tree & precalculate_distances) {
            
            for (i in seq_along(obs)) {
              
              if (!is.null(obs[[i]]$tree)) {
                obs_height = max(ape::cophenetic.phylo(obs[[i]]$tree)["GL",])
                mrate = target_height / obs_height
                lambda_edge = obs[[i]]$tree$edge.length * mrate
                
                set.seed(obs[[i]]$samples$seed)
                if ("dispersion" %in% rownames(simulation$params$param.matrix)) {
                  # sample edge length from negative-binomial distribution
                  disp = simulation$params$param.matrix["dispersion",1]
                  n_gens = obs[[i]]$tree$edge.length
                  mu = mrate
                  el = rnbinom(length(n_gens), mu=mu * n_gens, size=1/disp)
                } else {
                  # sample edge length from poisson (sum of poisson is poisson)
                  el = rpois(length(lambda_edge), lambda_edge)
                }
                obs[[i]]$tree$edge.length = el
                
              } else {
                mrate = NA
              }
              
              obs[[i]]$mutation_rate = mrate
              
            }
            
            mean_mrate = mean(sapply(obs, "[[", "mutation_rate"), na.rm = TRUE)
            particle$mutation_rate = mean_mrate
          }
          
        }
        
        
        # store all observations to a file
        if (!is.null(dump_dir)) {
          tryCatch({
            if (!is.null(generator_sample_params)) {
              
              stored_particle = particle
              if (is.null(simulation)) {
                stored_particle$observations = obs[1]
              } else {
                stored_particle$observations = obs
              }
              
              wh_hash = !names(stored_particle) %in% c("observations","verbose")
              hash = digest::digest(stored_particle[wh_hash])
              timestamp = as.hexmode(round(as.numeric(Sys.time())))
              out_file = file.path(dump_dir, paste0(hash, ".", timestamp, ".rds"))
              
              saveRDS(stored_particle, out_file, version = 2) 
            }
          }, error=function(e) return(NULL))
        }
        
        
        # calculate distances:
        if (precalculate_distances) {
          dists = sapply(obs, rho_, y=target_tree)
          if (!keep_trees) {
            for (i in seq_along(obs)) {
              obs[[i]]$tree = NULL
            }
          }
        }
        
        # only keep the oberservation closest to the target
        if (keep_min_obs_sim) {
          wh_min = which.min(dists)
          if (length(wh_min) == 0) wh_min = 1
          obs = obs[wh_min[1]]
          dists = dists[wh_min[1]]
        }
        
        # store pre-calculated distances in particle
        if (precalculate_distances) {
          particle$distances[[k]] = dists
        }
        
        
        if (length(obs) > 1) {
          particle$observations[[k]] = obs
        } else {
          particle$observations[[k]] = obs[[1]]
        }
        
      }
      
      return(particle)
    }, mc.cores = n_cores, mc.preschedule=FALSE)
    
    options(width=old_width) # revert print width
    
    return(result)
  }

sample_particles = 
  function(constant_params, generator_sim_params, N) {
    
    # fixed params for simulations:
    
    lapply(seq_len(N), function(i) {
      # generate particle
      particle = c(constant_params, generator_sim_params()) # random particle
      particle = particle[!duplicated(names(particle))]
      return(particle)
      
    })
    
  }




tree_distance = function(a, b, method="treeVec", swap_labels=FALSE, label_group_function=NULL, return_tree=FALSE, ...) {
  
  if (method %in% c("sEPD","EPD")) {
    .dist_mat = 
      function(cp1, cp2) {
        cp2 = cp2[rownames(cp1), colnames(cp1)]
        sqrt(sum((cp1[upper.tri(cp1)]-cp2[upper.tri(cp2)])^2))
      }
  } else if (method %in% c("sMPD","MPD")) {
    .dist_mat = 
      function(cp1, cp2) {
        cp2 = cp2[rownames(cp1), colnames(cp1)]
        sum(abs(cp1[upper.tri(cp1)]-cp2[upper.tri(cp2)]))
      }
  } else if (method %in% c("sCPD","CPD")) {
    .dist_mat = 
      function(cp1, cp2) {
        cp2 = cp2[rownames(cp1), colnames(cp1)]
        d_pairs = rbind(c(cp1[upper.tri(cp1)]),c(cp2[upper.tri(cp2)]))
        1-cor(d_pairs[1,],d_pairs[2,])
      }
  } else {
    stop("Invalid distance.")
  }
  
  .swap_label = function(d, i, j) {
    x = rownames(d)[i]
    colnames(d)[i] = rownames(d)[i] = rownames(d)[j]
    colnames(d)[j] = rownames(d)[j] = x
    return(d)
  }
  
  
  cp1 = ape::cophenetic.phylo(a)
  cp2 = ape::cophenetic.phylo(b)
  
  if (method %in% c("sEPD", "sMPD", "sCPD")) {
    # scale
    cp1 = cp1 / max(cp1[, "GL"], na.rm = TRUE) * 1000
    cp2 = cp2 / max(cp2[, "GL"], na.rm = TRUE) * 1000
  }
  
  cp2 = cp2[rownames(cp1), colnames(cp1)]
  c_dist = .dist_mat(cp1, cp2) # current distance
    
    
    if (swap_labels) {
      
      # indentify group of labels that can be swapped
      labels = rownames(cp1)
      
      if (is.null(label_group_function)) {
        label_groups = rep(1, length(labels))
      } else {
        label_groups = label_group_function(labels)
      }
      label_groups[labels == "GL"] = NA
      
      swappable = split(labels, label_groups)

      potential_swaps = 
        swappable %>% 
        lapply(function(x) expand.grid(x, x)) %>% 
        do.call(what=rbind) %>% 
        dplyr::filter(Var1 != Var2) %>% 
        apply(1, function(x) sort(x)) %>% t() %>% unique() %>% # drop duplicated swaps (i.e A->B, B->A)
        t() %>% data.frame() %>% as.list %>% 
        lapply(function(x) match(x, rownames(cp1)))
        

      # test actual swaps:
      if (length(potential_swaps)) {
          
        while (TRUE) {
          # swapped distances
          alt_cp2 = lapply(potential_swaps, function(x) .swap_label(cp2, x[1], x[2]))
          alt_dist = sapply(alt_cp2, .dist_mat, cp1=cp1)
          delta_d = alt_dist - c_dist
          
          wh_swap = which.min(delta_d)
          if (delta_d[wh_swap] >= 0) {
            break
          }
          
          cp2 = alt_cp2[[wh_swap]]
          c_dist = alt_dist[wh_swap]
        }
      }
    }
   
  if (return_tree) {
    swap_table = magrittr::set_names(colnames(cp2), colnames(cp1))
    b$tip.label = as.character(swap_table[b$tip.label])
    return(list(tree=b, dist=as.numeric(c_dist)))
  }
  
  return(as.numeric(c_dist))
}


continue_chain = function(dir, ...) {
  
  # find last state file
  particle_files = ordered_particle_list(dir)
  last_state_file = tail(particle_files, 1)
  
  # catch wrong calls
  stopifnot(file.exists(file.path(dir, "smc_chain", "options.rds")))
  
  # number of existing chains
  n_chain = length(list.files(dir, "smc_chain")) 
  out_dir_new_chain = file.path(dir, paste0("smc_chain", n_chain + 1))
  stopifnot(!file.exists(out_dir_new_chain))
  
  # load and modify options
  options = readRDS(file.path(dir, "smc_chain", "options.rds"))
  options$output_dir = out_dir_new_chain # new output dir.
  
  # optionally modify parameters:
  new_args = list(...)
  if (length(new_args)) {
    for (el in names(new_args)) {
      options[[el]] = new_args[[el]]
    }
  }
  
  # revert last state if existing:
  if (length(last_state_file) > 0) {
    
    vars_previous_state =
      c("particles" = "particles",
        "n_step" = "n_step",
        "acceptance_rate" = "avg_acceptance",
        "delta_d" = "delta_d",
        "epsilon" = "epsilon",
        "rng_state" = "rng_state"
      )
    
    last_state = readRDS(last_state_file)
    options$previous_state = last_state[names(vars_previous_state)]
    names(options$previous_state) = vars_previous_state
    
    # fill in data for number of steps and delta d if they are missing:
    if (is.null(options$previous_state$n_step)) {
      options$previous_state$n_step = sum(!grepl("step_0", particle_files))
    }
    
    if (is.null(options$previous_state$delta_d)) {
      wh_not_resume = !grepl("step_0", particle_files)
      eps_states = sapply(particle_files[wh_not_resume], function(x) readRDS(x)$epsilon)
      eps_states = c(Inf, eps_states)
      diff = (eps_states[-length(eps_states)] - eps_states[-1])
      rel_diff = as.numeric(diff / eps_states[-length(eps_states)])
      options$previous_state$delta_d = c(Inf, rel_diff)
    }
    
  }
  
  # call smc abc again
  smc_chain = do.call(smc_abc_tree, options)
  
  return(smc_chain)
}

