format_model_string = function(x, short=TRUE) {
  x = as.character(x)

  if (length(x) > 1)
    return(sapply(x, format_model_string))

  x = gsub("no_death_", "", x) %>%
    gsub(pattern = "selection_(?=2)", replacement = "selection x ", perl = TRUE) %>%
    gsub(pattern = "dispersion", replacement = "D") %>%
    gsub(pattern = "regionsize", replacement = "RS") %>%
    strsplit("_") %>%
    unlist() %>% Hmisc::capitalize() %>%
    paste0(collapse = " + ")

  if (short)  {
    x = gsub(" [+] Push", "", x)
    x = gsub(" [+] D( |$)", " ", x)
    x = gsub(" [+] RS( |$)", " ", x)
    x = gsub(" $", "", x)
  }

  x
}

smc_abc_model_selection = function(models, fname="result_set.rds", suffix="", title="", n_cores=1, ...) {

  ggsave = function(filename, plot, ...) {
    of = file.path(dirname(filename), paste0(".", gsub("[.][a-zA-Z]+$", ".rds", basename(filename))))
    saveRDS(plot, of)
    ggplot2::ggsave(filename, plot, ...)
  }
  
  ics_models = list()

  # define a resonable kernel for the surogate likelihood:
  in_files = file.path(models, fname)
  stopifnot(all(file.exists(in_files)))
  eps_min_smc = min(sapply(in_files, function(x) min(readRDS(x)$states$epsilon)))

  d = readRDS(in_files[1])
  eps_noise = get_average_epsilon_noisy_tree(d$target, d$variables$rho, level = 1)
  kernel_surrogate_lik = (function(d) { force(d); return(function(x) { y = x < d; storage.mode(y) = "numeric"; y }) })(eps_noise)

  for (i in seq_along(models)) {

    model = names(models)[i]
    cat("Evaluating model ", i, "/", length(models), " '", model, "':\n\n", sep="")

    tryCatch({

      # files
      in_file = file.path(models[model], fname)
      out_file = file.path(models[model], paste0("model_selection_data", suffix, ".rds"))
      fig_file = file.path(models[model], paste0("model_selection_data", suffix, ".pdf"))
      stopifnot(file.exists(in_file))

      # check if model file needs updating
      if (file.exists(out_file)) {
        if (file.mtime(out_file) < file.mtime(in_file)) {
          previous_model_selection_data = list()
        } else {
          previous_model_selection_data =
            tryCatch(readRDS(out_file), error=function(e) return(list()))
        }
      } else {
        previous_model_selection_data = list()
      }

      # update output if needed:
      new_results =
        get_aic_and_etc(
          dir = models[model],
          kernel_lik = kernel_surrogate_lik,
          results = previous_model_selection_data,
          fname=fname,
          n_cores = n_cores,
          ...
        )

      # add distance values
      new_results$eps = eps_noise
      new_results$eps_min_smc = eps_min_smc
      new_results$eps_noise = eps_noise

      # add synethic likelihood fits:
      cat("Adding syn-lik fits\n")
      new_results = add_synthetic_lik_fits(new_results, n_cores = n_cores)
      stats_sl = as.list(calc_model_stats_syn_lik(new_results, new_results$eps))
      names(stats_sl) = paste0(names(stats_sl), "_sl")
      new_results[names(new_results) %in% names(stats_sl) ] = NULL
      new_results = c(new_results, stats_sl)

      # export data (drop kernel due to large size of stored data)
      new_results_export = new_results
      new_results_export$kernel = NULL
      saveRDS(new_results_export, out_file)


      # plot 1/4) of distance under map
      try({
        plot =
          plot_model_dist_result(new_results, paste0(title, " (", format_model_string(model), ")")) +
          theme(plot.subtitle = element_text(size=9)) +
          theme(plot.caption = element_text(size=9)) +
          scale_x_continuous(n.breaks = 4)
        ggsave(fig_file, plot, width = 5.5, height = 3.5)
      })


      # plot 2/4) p value of model fit
      try({
        plot =
          plot_model_fit_estimate(new_results, paste0(title, " (", format_model_string(model), ")")) +
          scale_x_continuous(n.breaks = 4) +
          theme(plot.title = element_text(size=11))
        ofile = file.path(models[model], paste0("model_fit", suffix, ".pdf"))
        ggsave(ofile, plot, width = 3, height = 3.1)
      })


      # plot 3/4) qq plot of data
      try({
        set.seed(123)

        plot =
          normal_qq_plot_model_data(new_results, n_max = 50) +
          ggtitle(paste0(title, " (", format_model_string(model), ")")) +
          theme(plot.title = element_text(size=12))
        ofile = file.path(models[model], paste0("model_qqplot", suffix, ".pdf"))
        ggsave(ofile, plot, width=3.6, height=2.5)
      })

      # plot 4/4) plot of sl fits
      try({
        plot =
          plot_model_dist_result_sl(
            new_results,
            eps = new_results$eps,
            allow_sl = TRUE,
            paste0(title, " (", format_model_string(model), ")"),
            n_cores = n_cores
          ) + scale_x_continuous(n.breaks = 4)
        ofile = file.path(models[model], paste0("model_selection_data", suffix, "_sl.pdf"))
        ggsave(ofile, plot, width = 6.5, height = 3)
      })


      # drop some large elements and append to model selection dataset
      drop_el = grep("^d_", names(new_results), value=TRUE)
      drop_el = c(drop_el, grep("^particle", names(new_results), value=TRUE))
      drop_el = c(drop_el, grep("min_sc_frac$", names(new_results), value=TRUE))
      drop_el = c(drop_el, grep("trees$", names(new_results), value=TRUE))
      new_results$kernel = NULL
      ics_models[[model]] = new_results[!names(new_results) %in% drop_el]

      # print current results
      wh_print = c("AIC", "AIC2", "neg_log_lik","n_param","map","median","mean")
      print(ics_models[[model]][wh_print])
    }, error=function(e) print(e))

  }

  return(ics_models)
}

get_average_epsilon_noisy_tree = function(t, dist, level=1) {
  t_dash =
    lapply(1:1000, function(i) {
      noise = 2^runif(length(t$edge.length), -log2(1+level), log2(1+level))
      t$edge.length = t$edge.length * noise
      return(t)
    })

  dists = sapply(t_dash, function(x) dist(x, t))
  mean(dists)
}

get_aic_and_etc = function(dir, kernel_lik, fname, eps=NULL, N_p=0, N=0, N_sim=0, M=0, results=list(), n_cores=1, keep_trees=TRUE) {

  stopifnot(is.list(results))
  old_width = getOption("width")
  options(width=50)

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Wipe previous data if needed  ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (!all(c("N_p","N", "N_sim", "M") %in% names(results))) {
    results = list()
  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Load data ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  print_header = function(x) cat(crayon::bold(crayon::black(crayon::bgWhite(x))))
  print_subheader = function(x) cat(crayon::underline(x))

  cat("Loading dataset ... ")
  dump_file = file.path(dir, fname)
  result_set = readRDS(dump_file)

  # reassign weights:
  if (is.null(eps)) eps = min(result_set$states$epsilon, na.rm=TRUE)
  result_set = recalculate_fracs_below_eps(result_set, eps)
  sum_fracs = sum(result_set$particles$frac_below_eps, na.rm=TRUE)
  result_set$particles$weight = result_set$particles$frac_below_eps / sum_fracs
  result_set$particles$weight[is.na(result_set$particles$weight)] = 0

  cat(crayon::green("Done.\n\n"))

  # variables from result set
  target = MLLPT::set_lp_tiplength(result_set$variables$target_y, 0)
  rho = result_set$variables$rho
  s_gen = result_set$variables$s_gen
  tree_mod = result_set$variables$tree_manipulator

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Calculate MAP ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  print_header("-=-=-=-= Calculating MAP -=-=-=-=\n\n")
  map_data = get_map_particle(result_set, verbose = TRUE)

  # assume input data changed if the mean is different
  mean_n = results[["mean"]][!names(results[["mean"]]) %in% "mutation_rate"]
  mean_o = map_data[["mean"]][!names(map_data[["mean"]]) %in% "mutation_rate"]
  if ((!isTRUE(all.equal(mean_n, mean_o)))) {
    cat("  =>  MAP estimate changed! Updating all data!\n")
    results = list()
  }


  # update data if point estimate changed
  dhat_n = results[["est_dhat"]][!names(results[["est_dhat"]]) %in% "mutation_rate"]
  dhat_o = map_data[["est_dhat"]][!names(map_data[["est_dhat"]]) %in% "mutation_rate"]
  if ((!isTRUE(all.equal(dhat_n, dhat_o)))) {
    cat("  =>  Dhat point estimate changed! Updating dhat data!\n")
    print(cbind(results[["est_dhat"]], map_data[["est_dhat"]]))
    results[["d_dhat"]] = list()
  }

  wh_copy = names(map_data)
  results[wh_copy] = map_data[wh_copy]
  particles_dhat = rep(list(map_data[["particle_dhat"]]), N_p)

  cat("\n")

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # DHAT ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  print_header("-=-=-=-= Estimation of dhat -=-=-=-=\n\n")

  dhat_needs_update =
    is.null(results[["d_dhat"]]) | # missing data
    !is.numeric(results[["d_dhat"]]) | # bad data
    is.null(results[["d_dhat_permutated"]]) | # missing data
    !isTRUE(results[["N_sim"]] >= N_sim) | # increased number of simulations
    !isTRUE(results[["N_p"]] >= N_p) # increased number of simulations

  if (dhat_needs_update) {

    print_subheader(" => 1/3) Sampling D* under MAP\n")

    set.seed(123)

    obs_map =
      observe_particles(
        particles_dhat,
        s_gen,
        M = 1,
        M_sim = N_sim,
        n_cores = n_cores,
        sample_tree = TRUE, #ifelse(is.null(result_set$variables$sample_tree), TRUE, result_set$variables$sample_tree),
        change_seed = TRUE,
        min_frac_clones = 0,
        tree_manipulator = tree_mod
      )

    cat("\n")
    print_subheader(" => 2/3) Calculating distances\n")

    results[["d_dhat"]] =
      d_matrix_from_particles_pb(
        obs_map,
        target,
        rho,
        n_cores = n_cores
      )

    results[["d_dhat_min_sc_frac"]] =
      lapply(obs_map, function(x)
        lapply(x$observations[[1]], function(y) {
          if (length(y$cell_counts) == 1) return(NA)
          sc = min(y$cell_counts[-1]/sum(y$cell_counts))
        }) %>% do.call(what=c)
      ) %>% do.call(what=rbind)

    if (keep_trees){
      results[["d_dhat_trees"]] =
        lapply(obs_map, function(x)
          lapply(x$observations[[1]], function(y) {
            y$tree
          })
        )
    }


    print_subheader(" => 3/3) Calculating expected distances\n")

    # also calculate distance between trees under MAP
    trees =
      lapply(obs_map, function(x)
        lapply(x$observations[[1]], function(y) {
          y$tree
        })
      )

    results[["d_dhat_permutated"]] =
      do.call(rbind, pbmcapply::pbmclapply(seq_along(trees), function(i) {
        do.call(rbind, lapply(seq_len(5), function(z) {
          x = trees[[i]][[sample(length(trees[[i]]), 1)]]
          y_opts = do.call(c, trees[-i])
          sapply(length(trees[[1]]), function(j) { rho(x, y_opts[[sample(length(y_opts), 1)]]) })
        }))
      }, mc.cores = n_cores))

    cat("\n")
    stopifnot(is.numeric(results[["d_dhat"]]))

  } else {
    print_subheader(" => Using previously calculated data\n")
  }

  dbar_needs_update =
    is.null(results[["d_dbar"]]) |
    !is.numeric(results[["d_dbar"]]) |
    !isTRUE(results[["N_sim"]] >= N_sim) |
    !isTRUE(results[["N"]] >= N) |
    !isTRUE(results[["M"]] >= M)

  if (dbar_needs_update) {

    print_subheader(" => 1/2) Sampling D* from posterior\n")

    set.seed(123)

    wh = sample(seq_len(NROW(result_set$observations)), ifelse(M == 0, 2, M), replace = TRUE)
    idx = as.list(data.frame(t(result_set$observations[wh,c("i","j","k","l")])))
    obs_exp_DBAR = lapply(idx, get_sim_params, res = result_set)
    obs_exp_DBAR = lapply(obs_exp_DBAR, function(x) {x$seed = NA; return(x)})
    results[["particles_dbar"]] = obs_exp_DBAR

    obs_exp_DBAR =
      observe_particles(
        obs_exp_DBAR,
        s_gen,
        M = N,
        M_sim = ifelse(M == 0, 2, N_sim),
        sample_tree = ifelse(is.null(result_set$variables$sample_tree), TRUE, result_set$variables$sample_tree),
        n_cores = n_cores,
        change_seed = TRUE,
        min_frac_clones = 0,
        tree_manipulator = tree_mod
      )


    print_subheader(" => 2/2) Calculating distances\n")
    results[["d_dbar"]] =
      d_matrix_from_particles_pb(
        obs_exp_DBAR,
        target,
        rho,
        n_cores=n_cores
      )

    results[["d_dbar_min_sc_frac"]] =
      lapply(obs_exp_DBAR, function(x)
        lapply(x$observations, function(o){
          lapply(o, function(y) {
            if (length(y$cell_counts) == 1) return(NA)
            sc = min(y$cell_counts[-1]/sum(y$cell_counts))
          }) %>% unlist()
        }) %>% unlist()
      ) %>% do.call(what=rbind)

    if (keep_trees) {
      results[["d_dbar_trees"]] =
        lapply(obs_exp_DBAR, function(x)
          lapply(x$observations, function(o){
            lapply(o, function(y) { y$tree })
          })
        )
    }

    cat("\n")
    stopifnot(is.numeric(results[["d_dbar"]]))

  } else {
    print_subheader(" => Using previously calculated data\n")
  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # DBAR2 ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  print_header("-=-=-=-= Estimation of DBAR2 -=-=-=-=\n\n")

  dbar2_needs_update =
    is.null(results[["d_dbar2"]]) |
    !is.numeric(results[["d_dbar2"]]) |
    !isTRUE(results[["N_sim"]] >= N_sim) |
    !isTRUE(results[["N_p"]] >= N_p)

  if (dbar2_needs_update) {

    print_subheader(" => 1/2) Sampling D* from posterior\n")

    set.seed(123)

    wh = sample(seq_len(NROW(result_set$observations)), N_p, replace = TRUE)
    idx = as.list(data.frame(t(result_set$observations[wh,c("i","j","k","l")])))
    particles_exp_DBAR2 = lapply(idx, get_sim_params, res = result_set)
    particles_exp_DBAR2 = lapply(particles_exp_DBAR2, function(x) {x$seed = NA; return(x)})
    results[["particles_dbar2"]] = particles_exp_DBAR2

    obs_exp_DBAR2 =
      observe_particles(
        particles_exp_DBAR2,
        s_gen,
        M = 1,
        M_sim = N_sim,
        sample_tree = ifelse(is.null(result_set$variables$sample_tree), TRUE, result_set$variables$sample_tree),
        n_cores = n_cores,
        change_seed = TRUE,
        min_frac_clones = 0,
        tree_manipulator = tree_mod
      )

    print_subheader(" => 2/2) Calculating distances\n")

    results$d_dbar2 =
      d_matrix_from_particles_pb(
        obs_exp_DBAR2,
        target,
        rho,
        n_cores=n_cores
      )

    results[["d_dbar2_min_sc_frac"]] =
      lapply(obs_exp_DBAR2, function(x)
        lapply(x$observations[[1]], function(y) {
          if (length(y$cell_counts) == 1) return(NA)
          sc = min(y$cell_counts[-1]/sum(y$cell_counts))
        }) %>% do.call(what=c)
      ) %>% do.call(what=rbind)


    if (keep_trees) {
      results[["d_dbar2_trees"]] =
        lapply(obs_exp_DBAR2, function(x)
          lapply(x$observations[[1]], function(y) { y$tree })
        )
    }


    stopifnot(is.numeric(results$d_dbar2))
    cat("\n")

  } else {
    print_subheader(" => Using previously calculated data\n")
  }

  cat("\n")

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Print results ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cat("Finished collecting data!\n")
  results = calc_model_selection_stats(results, kernel_lik)
  cat(" => AIC = ", results[["AIC"]], "\n")
  cat(" => AIC2 = ", results[["AIC2"]], "\n")

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Return set ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  options(width=old_width)

  results[c("M","N_p","N","N_sim","kernel")] =
    list(M, N_p, N, N_sim, kernel_lik)

  return(results)
}

get_map_particle = function(x, verbose=FALSE) {

  get_map_1d = function(x, weights=NULL) {
    set.seed(123)
    wh = !is.na(x)
    kde = density(x[wh], weights = weights[wh], n = 100)
    kde$x[which.max(kde$y)]
  }


  pp = get_varied_parameter_list(x)
  wh = !colnames(pp) %in% c("weight","particle","weight_i","i","frac_below_eps")
  w = pp$weight
  w[is.na(w)] = 0
  pp = pp[,wh,drop=FALSE]
  result = list()

  stopifnot(is.numeric(w))

  # kde across all dims
  wh_tied = sapply(seq_len(NCOL(pp)), function(i) which(sapply(seq_len(NCOL(pp)), function(j) all(pp[,i] == pp[,j], na.rm=TRUE)))[1])
  pp_not_dup = pp[,unique(wh_tied),drop=FALSE]
  pp_not_dup = pp_not_dup[,colnames(pp_not_dup) != "mutation_rate", drop=FALSE]
  wh_use = !apply(is.na(pp_not_dup), 1, any)


  if (NCOL(pp_not_dup) < 5) {

    kde_all =
      ks::kde(
        pp_not_dup[wh_use,],
        w = w[wh_use] / sum(w) * sum(wh_use),
        xmin = apply(pp_not_dup, 2, min, na.rm = TRUE),
        xmax = apply(pp_not_dup, 2, max, na.rm = TRUE),
        gridsize = 100,
        binned = NCOL(pp_not_dup) < 5
      )


    wh_max_d = which(kde_all$estimate == max(kde_all$estimate), arr.ind=TRUE)

    if (is.null(nrow(wh_max_d))) {
      idx = as.numeric(wh_max_d[1])
    } else {
      idx = as.numeric(wh_max_d[1,])
    }

    if (length(idx) > 1) {
      map_kde = sapply(seq_along(idx), function(i) kde_all$eval.points[[i]][idx[i]])
    } else {
      map_kde = kde_all$eval.points[[idx]]
    }

    names(map_kde) = colnames(pp_not_dup)
    map_kde = map_kde[colnames(pp)[wh_tied]]
    names(map_kde) = colnames(pp)

  } else {
    map_kde = apply(pp, 2, get_map_1d, weights=w)
  }

  # get summary statistics:
  result$map_ndim = map_kde
  result$map = apply(pp, 2, get_map_1d, weights=w)
  result$map_ndim["mutation_rate"] = result$map["mutation_rate"]
  result$mean = apply(pp, 2, Hmisc::wtd.mean, weights=w,  normwt=TRUE, na.rm=TRUE)
  result$median = apply(pp, 2, Hmisc::wtd.quantile, weights=w, probs=0.5, normwt=TRUE, na.rm=TRUE)
  result$lower = apply(pp, 2, Hmisc::wtd.quantile, weights=w, probs=0.05, normwt=TRUE, na.rm=TRUE)
  result$upper = apply(pp, 2, Hmisc::wtd.quantile, weights=w, probs=0.95, normwt=TRUE, na.rm=TRUE)

  #set point estimates for particles
  result$est_dhat =  result$map_ndim

  # print values:
  if (verbose) {
    cat("Estimates:\n")

    estimate_string = paste0(
      names(result$map_ndim), ": ",
      signif(result$map_ndim, 3), " (",
      signif(result$lower, 3),"-",
      signif(result$upper, 3),")"
    )

    for (i in seq_along(estimate_string)) {
      cat("  -", estimate_string[i], "\n")
    }
  }

  # generate particles from map
  .get_particle_with_params = function(params) {
    params = params[!names(params) %in% "mutation_rate"] # ignore mutation rate
    idx = as.numeric(x$observations[1,c("i","j","k","l")])
    particle = get_sim_params(x, idx)
    for (i in seq_along(params)) {
      idx_row = gsub("[.].*", "", names(params[i]))
      idx_col = as.numeric(gsub(".*[.]", "", names(params[i])))
      particle$param.matrix[idx_row, idx_col] = params[i]
    }
    particle$observations = NULL
    particle$frac_below_eps = NULL
    return(particle)
  }

  result$particle_dhat = .get_particle_with_params(result$map_ndim)

  return(result)
}

plot_model_selection_result = function(msd, title="Model selection", stats=c("n_param","neg_log_lik","AIC", "AIC2"), alt_labels=c()) {

  # alternative labels:
  alt_labels_default = c(
    "n_param" = "k",
    "neg_log_lik" = "NLL",
    "neg_log_lik2" = "NLL'",
    "AIC2" = "AIC'"
  )

  wh_default = !names(alt_labels_default) %in% names(alt_labels)
  alt_labels = c(alt_labels, alt_labels_default[wh_default])

  stats_alt = stats
  wh = stats_alt %in% names(alt_labels)
  stats_alt[wh] = alt_labels[stats_alt[wh]]


  # get current subset of data:
  k = sapply(msd, "[[", "n_param")
  eps_vals = sapply(msd, function(x) ifelse("eps" %in% names(x), x$eps, NA))
  msd = lapply(msd, function(x) x[names(x) %in% stats])
  mean_vals = lapply(msd, lapply, mean, na.rm=TRUE)

  # find stats with bootstrap values:
  has_bs_all = lapply(lapply(msd, sapply, length), ">", 1)
  has_bs = has_bs_all[[1]]
  for (i in seq_along(has_bs_all))
    stopifnot(isTRUE(all.equal(has_bs_all[[i]][names(has_bs)], has_bs)))


  # restructure data
  model_selection_rs =
    reshape::melt(mean_vals) %>%
    dplyr::filter(!L2 %in% names(which(has_bs))) %>%
    dplyr::mutate(L2 = factor(L2, levels = stats, labels = stats_alt))


  model_selection_bs =
    reshape::melt(msd) %>%
    dplyr::filter(L2 %in% names(which(has_bs))) %>%
    dplyr::mutate(L2 = factor(L2, levels = stats, labels = stats_alt))


  # identify best values
  min_value_ic = with(reshape2::melt(mean_vals), tapply(value, L2, min))
  model_selection_rs$best = with(model_selection_rs, value == min_value_ic[as.character(L2)])
  model_selection_rs$best[model_selection_rs$L2 %in% c("k","n_param")] = FALSE
  model_selection_rs$best = FALSE

  if (NROW(model_selection_bs)) {
    model_selection_bs$best = with(model_selection_bs, value == min_value_ic[as.character(L2)])
    model_selection_bs$best[model_selection_rs$L2 %in% c("k","n_param")] = FALSE
    model_selection_bs$best = FALSE
  }

  # order model by number of parameters
  order_models = names(sort(k))
  new_level_model = format_model_string(order_models)
  model_selection_rs$L1 = factor(model_selection_rs$L1, order_models, new_level_model, ordered = TRUE)
  model_selection_bs$L1 = factor(model_selection_bs$L1, order_models, new_level_model, ordered = TRUE)

  # line for some:
  data_hlines =
    model_selection_rs %>%
    dplyr::filter(L2 %in% stats_alt) %>%
    dplyr::filter(L1 == "Neutral") %>%
    dplyr::mutate(value=value-5) %>%
    dplyr::filter(value >= 0)

  # create plot of results
  plot_model_selection =
    model_selection_rs %>%
    ggplot(aes(x=L1, y=value, fill=best)) +
    facet_wrap(~L2, nrow=1, scales = "free_y") +
    geom_bar(stat="identity", width=0.85) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    xlab("") +
    ylab("Value") +
    guides(fill="none") +
    scale_fill_manual(values=c("FALSE"="gray40", "TRUE"="#C21807")) +
    #geom_hline(data=data_hlines, aes(yintercept=value), linetype=3) +
    ggtitle(title)

    if (NROW(model_selection_bs)) {
      plot_model_selection =
        plot_model_selection +
        geom_violin(data=model_selection_bs, scale="width") +
        stat_summary(data=model_selection_bs, fun=mean)
    }

  if (!all(is.na(eps_vals))) {
    eps_vals = unique(na.omit(eps_vals))
    if (length(eps_vals) == 1) {
      plot_model_selection =
        plot_model_selection +
        labs(captions=bquote(epsilon==.(signif(eps_vals, 3))))
    }
  }

  return(plot_model_selection)
}

add_bootstrap_ms_statistics = function(x, eps=NA, n=100, n_cores=1) {

  stats_to_bootstrap = c("Dbar","Dbar2","Dhat", "AIC","AIC2",
                         "neg_log_lik2","neg_log_lik")

  stats_to_bootstrap = stats_to_bootstrap[stats_to_bootstrap %in% names(x)]

  if (is.na(eps)) {
    stopifnot("eps" %in% names(x))
    eps = x[["eps"]]
  }

  kernel = (function(d) {
    force(d)
    return(function(x) {
      y = x < d
      storage.mode(y) = "numeric"
      y
    })
  })(eps)

  bs_stats =
    pbmcapply::pbmclapply(seq_len(n), function(i) {
      unlist(.bootstrap_stats(x, calc_model_selection_stats, kernel_lik=kernel)[stats_to_bootstrap])
    }, mc.cores = n_cores) %>% do.call(what=rbind)

  for (stat in stats_to_bootstrap) {
    x[[stat]] = bs_stats[,stat]
  }

  return(x)
}

add_bootstrap_ms_statistics_syn = function(x, eps=NA, n=100, n_cores=1) {

  stats_to_bootstrap = c("Dbar","Dbar2","Dhat","AIC","AIC2", "neg_log_lik2","neg_log_lik")

  stats_to_bootstrap = stats_to_bootstrap[stats_to_bootstrap %in% names(x)]
  if (is.na(eps)) {
    stopifnot("eps" %in% names(x))
    eps = x[["eps"]]
  }

  cat("Bootstraping syn. lik. data.\n")
  bs_stats =
    pbmcapply::pbmclapply(seq_len(n), function(i) {
      unlist(.bootstrap_stats(x, calc_model_stats_syn_lik, eps=eps)[,stats_to_bootstrap])
    }, mc.cores = n_cores) %>% do.call(what=rbind)

  for (stat in stats_to_bootstrap) {
    x[[stat]] = bs_stats[,stat]
  }

  return(x)
}

get_neutral_cross_obs_distance = function(dfile, n_cores=1, ...) {

  N_kernel=10
  N_sim_kernel=100

  cat("Estimating distances for kernel ...\n")
  res = readRDS(dfile)
  rho = res$variables$rho

  # neutral_params:
  neutral_params =
    structure(c(1, 0.5, 1, 1, 50, 0, 0, -1), .Dim = c(8L, 1L), .Dimnames = list(
      c("birthrates", "deathrates", "aggressions", "push_power",
        "mutation_rates", "clone_start_times", "kill_regrow_times",
        "father"), "1"))

  # generate particles from map
  idx = as.numeric(res$observations[1,c("i","j","k","l")])
  particle = get_sim_params(res, idx)
  particle$param.matrix = neutral_params
  particles = rep(list(particle), N_kernel)


  # observe trees:
  observations =
    observe_particles(
      particles,
      res$variables$s_gen,
      M = 1,
      M_sim = N_sim_kernel,
      n_cores = n_cores,
      sample_tree = res$variables$sample_tree,
      change_seed = TRUE,
      tree_manipulator = res$variables$tree_manipulator
    )


  # also calculate distance between trees under MAP
  trees =
    lapply(observations, function(x) {
      lapply(x$observations[[1]], function(y) {  y$tree })
    }) %>% do.call(what=c)

  permutated =
    pbmcapply::pbmclapply(seq_len(N_kernel*100), function(x) {
      i = sample(length(trees), 1)
      j = sample(length(trees), 1)
      tryCatch(rho(trees[[i]], trees[[j]]), error=function(e) return(NA))
    })

  return(unlist(permutated))
}

alt_sgen = function (simulation, args, seed_samples = NULL, phi_A=NULL) {

  for (i in seq_along(args)) assign(names(args)[i], args[[i]])


  if ("region_scale_factor" %in% rownames(simulation$params$param.matrix)) {
    sf = simulation$params$param.matrix["region_scale_factor", 1]
    region_diameter = ceiling(region_diameter * sf)
  }

  max_samples_per_region = max(table(sapply(sample_annot,  "[[", "region")))
  if (prod(region_diameter) < max_samples_per_region) stop("Sample regions too small.\n")
  if (is.null(seed_samples)) seed_samples = sample(.Machine$integer.max, 1)
  set.seed(seed_samples)

  if (is.null(phi_A)) phi_A = runif(1, 0, 2 * pi)
  center = with(simulation$params, c(x, y, z)/2)
  max_r = with(simulation$params, max(c(x, y, z)/2))
  r_test = seq(0, max_r, by = 0.25)
  if (is.null(alternative_offset_regions)) {
    offsets_regions = c(A = 0, B = 0.5, C = 1, D = 1.5)
  } else {
    offsets_regions = alternative_offset_regions
  }
  region_center = lapply(offsets_regions * pi, function(phi_o) {
    phi = phi_A + phi_o
    test_points = data.frame(x = floor(r_test * cos(phi) + center[1]) - 1,
                             y = floor(r_test * sin(phi) + center[2]) - 1,
                             r_test = r_test) %>% unique()
    for (i in rev(seq_len(NROW(test_points)))) {
      if (simulation$sim$CellType(test_points$x[i], test_points$y[i],
                                  0)) {
        break
      }
    }
    r = test_points$r_test[i] * edge_distance
    x = floor(r * cos(phi) + center[1])
    y = floor(r * sin(phi) + center[2])
    z = 1
    c(x, y, z)
  })

  .get_pos_vec = function(x) seq(floor(-(x - 1)/2), floor((x -  1)/2))
  pos_locs = do.call(expand.grid, lapply(region_diameter, .get_pos_vec))
  sample_regions = sapply(sample_annot, "[[", "region")
  samples_per_region = split(names(sample_annot), sample_regions)
  sample_centers = list()

  for (i in seq_along(samples_per_region)) {
    c_center = region_center[[names(samples_per_region)[i]]]
    n_needed = length(samples_per_region[[i]])
    pos_sampled = list()
    idx_locs_p = seq_len(NROW(pos_locs))
    while (length(pos_sampled) < n_needed) {
      idxs_proposed = sample(idx_locs_p, size = n_needed, replace = FALSE)
      idx_locs_p = idx_locs_p[!idx_locs_p %in% idxs_proposed]
      offset = pos_locs[idxs_proposed, , drop = FALSE]
      centers_proposed = split(t(t(offset) + c_center),  seq_len(NROW(offset)))
      .contains_cell = function(x) do.call(simulation$sim$CellType, as.list(x - 1)) != 0
      is_accepted = sapply(centers_proposed, .contains_cell)
      pos_sampled = c(pos_sampled, centers_proposed[is_accepted])
    }
    pos_sampled = pos_sampled %>% magrittr::set_names(samples_per_region[[i]])
    sample_centers = c(sample_centers, pos_sampled)
  }

  centers = sample_centers %>% do.call(what = rbind)
  sampling_setup = sample_annot
  for (i in seq_along(sampling_setup)) {
    seed = sample.int(.Machine$integer.max, 1)
    sampling_setup[[i]]$center = as.numeric(centers[names(sampling_setup)[i], ])
    sampling_setup[[i]] = c(sampling_setup[[i]], min_vaf = 0, seed = seed)
  }

  return(list(setup = sampling_setup, vars = list(phi = phi_A),
              seed = seed_samples, region_diameter = region_diameter,
              region_center = region_center))
}

get_neutral_cross_obs_distance2 = function(dfile, n_cores=1, ...) {

  N_sims = 20
  N_samples_per_sim = 25

  cat("Estimating distances for kernel ...\n")
  res = readRDS(dfile)
  rho = res$variables$rho

  # neutral_params:
  neutral_params =
    structure(c(1, 0.5, 1, 1, 50, 0, 0, -1), .Dim = c(8L, 1L), .Dimnames = list(
      c("birthrates", "deathrates", "aggressions", "push_power",
        "mutation_rates", "clone_start_times", "kill_regrow_times",
        "father"), "1"))

  # generate particles from map
  idx = as.numeric(res$observations[1,c("i","j","k","l")])
  particle = get_sim_params(res, idx)
  particle$param.matrix = neutral_params

  #
  dists = lapply(seq_len(N_sims), function(a) {

    sim = do.call(new_simulation, particle)
    vars = ls(envir = environment(res$variables$s_gen))
    vals = mget(vars, envir = environment(res$variables$s_gen))
    phi = runif(1, 0, 2*pi)

    trees =
      pbmcapply::pbmclapply(seq_len(N_samples_per_sim), function(i){
        attach(vals)
        sample_annot = vals$sample_annot
        samples = alt_sgen(sim, phi_A = phi)
        get_equivalent_single_cell_tree(sim$sim, samples$setup)
      }, mc.cores = n_cores)


    permutated =
      pbmcapply::pbmclapply(seq_along(trees), function(i) {
        sapply(seq_along(trees)[-i], function(j) {
          tryCatch(rho(trees[[i]], trees[[j]]), error=function(e) return(NA))
        })
      }, mc.cores = n_cores)

    unlist(permutated)
  })


  return(dists)
}

get_best_sim_cross_obs_distance = function(dfiles, n_cores=1, ...) {

  N_trees = 100

  cat("Estimating distances for kernel ...\n")
  res_list = lapply(dfiles, readRDS)
  dist_min = sapply(res_list, function(x) min(x$observations$dist))
  res = res_list[[which.min(dist_min)]]
  rho = res$variables$rho

  # generate particles from map
  idx = as.numeric(res$observations[1,c("i","j","k","l")])
  sim = get_sim_by_idx(res, idx, scale_tree = TRUE, add_samples = TRUE)
  vars = ls(envir = environment(res$variables$s_gen))
  vals = mget(vars, envir = environment(res$variables$s_gen))

  trees =
    pbmcapply::pbmclapply(seq_len(N_trees), function(i){
      samples = alt_sgen(sim, phi_A = sim$phi, args=vals)
      get_equivalent_single_cell_tree(sim$sim, samples$setup)
    }, mc.cores = n_cores)


  dists =
    pbmcapply::pbmclapply(seq_along(trees), function(i) {
      sapply(seq_along(trees)[-i], function(j) {
        tryCatch(rho(trees[[i]], trees[[j]]), error=function(e) return(NA))
      })
    }, mc.cores = n_cores)

  return(unlist(dists))
}

plot_model_dist_result = function(msd, title="") {

  dist_data =
    data.frame(
      d = as.numeric(c(msd$d_dhat_permutated, c(msd$d_dhat))),
      set = c(rep("Delta (D*'*\\'',D*'*')*','~D*'*'%~%p('.|'*bar(theta))", length(msd$d_dhat_permutated)),
              rep("Delta (D,D*'*')*','~D*'*'%~%p('.|'*bar(theta))", length(c(msd$d_dhat))))
    )# %>% dplyr::mutate(set=factor(set, c("MAP vs. MAP","D vs. MAP"), ordered = TRUE))


  dist_data_per_group = NULL
  for (i in seq_len(min(c(100, NROW(msd$d_dhat))))) {
    wh = !is.infinite(msd$d_dhat[i,])
    d = density(msd$d_dhat[i,wh], n = 100)
    dist_data_per_group =
      data.frame(
        i = i,
        d = d$x,
        dens = d$y / sum(d$y) * length(c(msd$d_dhat)),
        set = "Delta (D,D*'*')*','~D*'*'%~%p('.|'*bar(theta))"
      ) %>% rbind(what=dist_data_per_group)
  }

  for (i in seq_len(min(c(100, NROW(msd$d_dhat_permutated))))) {
    wh = !is.infinite(msd$d_dhat_permutated[i,])
    d = density(msd$d_dhat_permutated[i,wh], n = 100)
    dist_data_per_group =
      data.frame(
        i = i,
        d = d$x,
        dens = d$y / sum(d$y) * length(c(msd$d_dhat_permutated)),
        set = "Delta (D*'*\\'',D*'*')*','~D*'*'%~%p('.|'*bar(theta))"
      ) %>% rbind(what=dist_data_per_group)
  }


  if (!is.null(msd$kernel)) {
    wh = !is.infinite(dist_data$d)
    d_kernel = seq(0, max(dist_data$d[wh], na.rm=TRUE), length.out = 1000)
    d_kernel = data.frame(d = d_kernel, l = msd$kernel(d_kernel))
    d_kernel$l = d_kernel$l / sum(d_kernel$l) * NROW(dist_data) * 2
  } else {
    d_kernel = NULL
  }

  caption =
    paste0(
      "AIC = ", signif(mean(msd$AIC, na.rm=TRUE), 3),
      ", AIC' = ", signif(mean(msd$AIC2, na.rm=TRUE), 3)
    )

  subtitle =
    paste0(
      "p(D) = ", signif(exp(-msd$neg_log_lik), 3),
      ", k = ", msd$n_param
    )

  hist =
    dist_data %>%
    dplyr::mutate(set=factor(set, unique(dist_data$set))) %>%
    ggplot(aes(x=d, fill=set, color=set)) +
    geom_histogram(bins = 100, position="identity", alpha=0.8, color=NA) +
    labs(fill="Set", color="Set") +
    xlab("Distance") +
    ylab("Counts") +
    labs("") +
    geom_line(
      data=dist_data_per_group %>%
        dplyr::mutate(set=factor(set, unique(dist_data$set))),
      aes(y=dens, group=paste0(set,i)), alpha=0.1
    ) +
    scale_fill_brewer(palette="Set1", direction=-1, labels=scales::parse_format()) +
    scale_color_brewer(palette="Set1", direction=-1, labels=scales::parse_format()) +
    ggtitle(title, subtitle = subtitle) +
    scale_color_grey(labels=scales::parse_format()) +
    #scale_fill_grey() +
    labs(caption = caption)

  if (!is.null(d_kernel)) {
    hist = hist +
      geom_line(data=d_kernel, aes(x=d, y=l), color="black", inherit.aes=FALSE)
  }

  hist

  return(hist)
}

reduce_barcode = function(x, short=FALSE) {
  x = gsub("_[CDRL][0-9]$", "", gsub("EPICC_", "", x))
  if (short) x = gsub("^C[0-9]+_", "", x)
  return(x)
}

calc_exp_lik = function(x) {
  x[is.na(x) | is.infinite(x)] = 0
  if (all(x %in% c(0,1))) { # assume rejection kernel
    est = (sum(x) + 1e-12) / (length(x) + 1e-12)
  } else {
    est = mean(c(x), na.rm=TRUE)
  }
  return(est)
}

.bootstrap_stats = function(x, f, ...) {
  for (el in c("d_dbar","d_dbar2","d_dhat")) {
    if (el %in% names(x)) {
      idx = sample(seq_len(NROW(x[[el]])), replace = TRUE)
      x[[el]] = x[[el]][idx,]
    }
  }

  for (el in c("d_dbar_syn","d_dbar2_syn")) {
    if (el %in% names(x)) {
      idx = sample(seq_along(x[[el]]), replace = TRUE)
      x[[el]] = x[[el]][idx]
    }
  }

  f(x, ...)
}

calc_model_selection_stats = function(x, kernel_lik=NULL) {

  if (is.null(kernel_lik)) {
    kernel_lik = x$kernel
    stopifnot(!is.null(kernel_lik))
  }

  # AIC
  log_lik = log(calc_exp_lik(kernel_lik(x[["d_dbar"]])))
  x[["n_param"]] = length(unique(as.numeric(x[["mean"]])))
  x[["neg_log_lik"]] = -log_lik
  x[["AIC"]] = 2 * x[["n_param"]] - 2 * log_lik

  log_lik = log(calc_exp_lik(kernel_lik(x[["d_dbar2"]])))
  x[["n_param"]] = length(unique(as.numeric(x[["mean"]])))
  x[["neg_log_lik2"]] = -log_lik
  x[["AIC2"]] = 2 * x[["n_param"]] - 2 * log_lik

  p_exp_ = apply(kernel_lik(x[["d_dbar"]]), 1, calc_exp_lik)
  x[["Dbar"]] = mean(-2*log(p_exp_))
  x[["Dhat"]] =  -2 * log(calc_exp_lik(kernel_lik(x[["d_dhat"]])))

  p_exp_ = apply(kernel_lik(x[["d_dbar2"]]), 1, calc_exp_lik)
  x[["Dbar2"]] = mean(-2 * log(p_exp_))

  # p-value of model:
  mean_d = apply(x[["d_dhat"]] , 1, mean)
  mean_d_permut = apply(x[["d_dhat_permutated"]], 1, mean)
  x[["p_model"]] = mean(mean_d_permut > mean_d)

  return(x)
}

get_model_p_value = function(msd) {
  mean_d = apply(msd$d_dhat, 1, mean)
  mean_d_permut = apply(msd$d_dhat_permutated, 1, mean)
  signif(mean(mean_d_permut > mean(mean_d)), 3)
}

plot_model_fit_estimate = function(msd, title="") {

  # calculate mean dist
  mean_d = apply(msd$d_dhat, 1, mean)

  estimates_obs =
    data.frame(
      mean = mean(mean_d),
      low = quantile(mean_d, 0.01),
      high = quantile(mean_d, 0.99),
      set = c("MAP")
    )

  # calulate mean dist for permuated data:
  mean_d_permut = apply(msd$d_dhat_permutated, 1, mean)
  p_est = signif(mean(mean_d_permut > estimates_obs$mean), 3)
  subtitle = bquote(p(bar(Delta[sim]) > bar(Delta[obs]))%~~%.(p_est))

  hist =
    data.frame(x=mean_d_permut) %>%
    ggplot(aes(x=x)) +
    geom_histogram(bins = 100, position="identity", alpha=1, fill="gray10") +
    labs(fill="Set", color="Set") +
    xlab("Distance") +
    ylab("Counts") +
    labs("") +
    geom_vline(data=estimates_obs, aes(xintercept=mean), color="red", linetype=2) +
    ggtitle(title) +
    scale_color_grey(labels=scales::parse_format()) +
    labs(caption = subtitle)

  hist
}

add_synthetic_lik_fits = function(model_set, n_cores=n_cores) {

  .get_lik_fit =function(d) {
    d[is.infinite(d)] = NA
    if (sum(!is.na(d)) == 0) return(NA)
    fit = mclust::densityMclust(d, modelNames="V", verbose=FALSE)
    return(fit)
  }

  model_set$d_dhat_syn = .get_lik_fit(c(model_set$d_dhat))
  model_set$d_dbar_syn = pbmcapply::pbmclapply(seq_len(NROW(model_set$d_dbar)), function(i) .get_lik_fit(model_set$d_dbar[i,]), mc.cores = n_cores)
  model_set$d_dbar2_syn = pbmcapply::pbmclapply(seq_len(NROW(model_set$d_dbar2)), function(i) .get_lik_fit(model_set$d_dbar2[i,]), mc.cores = n_cores)

  return(model_set)
}

calc_model_stats_syn_lik = function(msd, eps, return_liks=FALSE, use_sl=TRUE, n_cores=1) {

  .get_lik =function(fit, eps) {

    if (isTRUE(all.equal(fit, NA))) return(NA)

    if ("Mclust" %in% class(fit)) {
      if (use_sl) {
        est = log(mclust::cdfMclust(fit, eps)$y)
      } else {
        est = mean(fit$data <= eps)
      }
    } else {
      est = mean(fit <= eps)
    }

    return(est)
  }

  s_lik_map = .get_lik(msd$d_dhat_syn, eps)
  s_lik_post = unlist(lapply(msd$d_dbar_syn, .get_lik, eps=eps))
  s_lik_post2 = unlist(lapply(msd$d_dbar2_syn, .get_lik, eps=eps))

  stats =
    data.frame(
      neg_log_lik = -log(mean(exp(s_lik_post), na.rm = TRUE)),
      neg_log_lik2 = -log(mean(exp(s_lik_post2), na.rm = TRUE)),
      neg_log_lik_map = -s_lik_map,
      k = length(msd$mean),
      Dbar = mean(-2 * s_lik_post, na.rm = TRUE),
      Dbar2 = mean(-2 * s_lik_post2,  na.rm = TRUE),
      Dhat =  -2 * s_lik_map,
      pV = var(-2 * s_lik_post, na.rm = TRUE) / 2,
      pV2 = var(-2 * s_lik_post2, na.rm = TRUE) / 2
    ) %>%
    dplyr::mutate(
      AIC = 2 * k - 2 * -neg_log_lik, 
      AIC2 = 2 * k - 2 * -neg_log_lik2
    )

  if (return_liks) {
    return(
      list(
        stats = stats,
        s_lik_map = s_lik_map,
        s_lik_post = s_lik_post,
        s_lik_post2 = s_lik_post2
      )
    )
  } else {
    return(stats)
  }
}

plot_model_dist_result_sl = function(msd, eps, allow_sl=TRUE, title="", n_cores=1) {

  # calculate density estimates:
  dist_data_per_group = NULL

  # a) map
  wh = !is.infinite(c(msd$d_dhat)) & !is.na(c(msd$d_dhat))
  dens = density(c(msd$d_dhat)[wh], n = 100)
  dist_data_per_group =
    data.frame(i = 1, d = dens$x, y = dens$y, set = "MAP(theta)") %>%
    rbind(what=dist_data_per_group)

  # b) posterior
  for (i in head(seq_len(NROW(msd$d_dbar2)), 200)) {
    wh = !is.infinite(msd$d_dbar2[i,]) & !is.na(msd$d_dbar2[i,])
    if (sum(wh) == 0) next()
    dens = density(msd$d_dbar2[i,wh], n = 100)
    dist_data_per_group =
      data.frame(i = i, d = dens$x, y = dens$y, set = "p(bar(theta)*'.|'*D)") %>%
      rbind(what=dist_data_per_group)
  }
  rownames(dist_data_per_group) = NULL

  hist =
    dist_data_per_group %>%
    ggplot(aes(x=d, y=y, color=set, group=paste0(i, set))) +
    cowplot::theme_cowplot() +
    geom_line(alpha=0.7) +
    geom_line(data=dist_data_per_group[dist_data_per_group$set=="MAP(theta)",], alpha=1) +
    labs(fill="", color="") +
    xlab("Distance") +
    ylab("Density") +
    scale_color_brewer(palette="Set1", direction=1, labels=scales::parse_format()) +
    ggtitle(title) +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(size=11)) +
    geom_vline(xintercept = eps, linetype=3)

  #
  if (!is.na(eps)) {


    syn_lik_data =
      calc_model_stats_syn_lik(
        msd,
        eps,
        return_liks = TRUE,
        use_sl = allow_sl,
        n_cores = n_cores
      )

    caption =
      with(syn_lik_data$stats, {
        paste0(
          "AIC = ", signif(mean(AIC, na.rm=TRUE), 3),
          ", AIC' = ", signif(mean(AIC2, na.rm=TRUE), 3)
        )
      })


    subtitle =
      with(syn_lik_data$stats, {
        paste0("p(D) = ", signif(exp(-neg_log_lik), 3), ", k = ", k)
      })


    hist2 =
      data.frame(x=syn_lik_data$s_lik_post) %>%
      ggplot(aes(x=x)) +
      cowplot::theme_cowplot() +
      geom_histogram(aes(fill="PP"), alpha=0.8, bins=100) +
      geom_histogram(data=data.frame(x=syn_lik_data$s_lik_post2), aes(fill="PP'"), alpha=0.75, bins=100) +
      geom_vline(data = syn_lik_data$stats, aes(xintercept=-neg_log_lik_map, linetype="MAP")) +
      xlab("Log-Likelihood") + ylab("Count") +
      labs(linetype="", fill="") +
      theme(legend.position = "bottom") +
      labs(caption=caption) +
      scale_fill_brewer(palette = "Set1", direction = -1)


    hist = cowplot::plot_grid(hist, hist2, align = "hv")
  }

  return(hist)
}

normal_qq_plot_model_data = function(msd, n_max=10) {

  qq_data = NULL

  idx = sample(seq_len(NROW(msd$d_dhat)), size = min(c(NROW(msd$d_dhat), n_max)), replace = FALSE)
  for (i in idx) {
    d = data.frame(qqnorm(msd$d_dhat[i,], plot.it = FALSE), set="f(z*'|'*hat(theta))", i=i)
    qq_data = rbind(qq_data, d)
  }

  idx = sample(seq_len(NROW(msd$d_dbar)), size = min(c(NROW(msd$d_dbar), n_max)), replace = FALSE)
  for (i in idx) {
    d = data.frame(qqnorm(msd$d_dbar[i,], plot.it = FALSE), set="f(z*'|'*theta) ~ p(theta*'|'*D)", i=i)
    qq_data = rbind(qq_data, d)
  }

  idx = sample(seq_len(NROW(msd$d_dbar2)), size = min(c(NROW(msd$d_dbar2), n_max)), replace = FALSE)
  for (i in idx) {
    d = data.frame(qqnorm(msd$d_dbar2[i,], plot.it = FALSE), set="f(z*'|'*theta) ~ p(theta*'|'*D)", i=i)
    qq_data = rbind(qq_data, d)
  }


  plot =
    ggplot(qq_data, aes(x=x, y=y, group=i)) +
    geom_line(alpha=0.7) +
    facet_wrap(~set, labeller = label_parsed) +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    ggtitle("Normal Q-Q Plot") +
    geom_smooth(aes(group=set), formula = y~x, method="lm", se=FALSE, color="red", linetype=2, size=0.7)

}
