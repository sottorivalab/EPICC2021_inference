load_smc_result_set = function(dir, include_rejected_particles=TRUE, n_cores=1) {
  
  print(dir)
  
  # helper function
  param_matrix_to_data.frame = 
    function(params) {
      params %>% 
        reshape::melt() %>% 
        dplyr::mutate(parameter=paste0(X1, ".", X2)) %>% 
        dplyr::select(-X1, -X2) %>% 
        reshape2::acast(.~parameter, value.var = "value") %>%
        as.data.frame()
    }
  
  load_particle_set = 
    function(x, eps=NA) {
      
      data = readRDS(x)
      
      if (all(c("epsilon","ess","acceptance_rate") %in% names(data))) {
        state_data = as.data.frame(data[c("epsilon","ess","acceptance_rate")])
        if (is.na(eps)) eps = data$epsilon
      } else { # format of rejected particles
        state_data = NULL
        data = list(particles=data)
        data$weight = rep(NA, length(data$particles))
        if (is.na(eps)) stop("Missing eps for rejected particles.")
      }
      
      dists = d_matrix_from_particles(data$particles, NULL, NULL, simplify=TRUE)     
      if (is.null(dim(dists))) return(NULL)
      new_weight = calc_w(dists, eps)
      new_weight[is.na(new_weight)] = 0
      
      # get particle data
      wh_cols = names(data$particles[[1]])
      wh_cols_not = c("param.matrix","observations","distances","seed")
      wh_cols = wh_cols[!wh_cols %in% wh_cols_not]
      
      particle_vars = 
        lapply(data$particles, "[", wh_cols) %>% 
        do.call(what=rbind) %>% 
        as.data.frame()
      
      for (i in seq_len(NCOL(particle_vars))) {
        particle_vars[,i] = 
          sapply(particle_vars[,i], function(x){if (is.null(x)) NA else x}) %>% 
          unlist()
      }
      
      particle_params = 
        lapply(data$particles, "[", "param.matrix") %>%
        lapply(param_matrix_to_data.frame) %>%
        do.call(what=rbind) %>%
        magrittr::set_rownames(NULL)
      
      if ("mutation_rate" %in% names(data$particles[[1]])) {
        particle_params$mutation_rate_sf = 
          sapply(data$particles, function(x) ifelse(is.null(x$mutation_rate), NA, x$mutation_rate))
      }
      
      
      particle_data = cbind(
        j = seq_along(data$weight),
        particle_vars,
        particle_params,
        weight_i = data$weight,
        weight = new_weight
      )
      
      
      # realisation data
      null_to_na = function(x) ifelse(is.null(x), NA, x)
      particle_obs_list = lapply(data$particles, "[[", "observations")
      first_obs = lapply(particle_obs_list, lapply, "[[", 1)
      
      sim_seeds = lapply(first_obs, lapply, "[[", "seed") 
      sim_seeds = lapply(sim_seeds, lapply, null_to_na) %>% reshape::melt()
      sim_seeds = sim_seeds[,c("L1","L2","value")]
      colnames(sim_seeds) = c("j","k","seed")
      
      cell_counts = lapply(first_obs, lapply, "[[", "cell_counts")
      cell_counts = lapply(cell_counts, lapply, function(x) { if(is.null(x)) return(NA); magrittr::set_colnames(data.frame(sum(x), matrix(x, nrow=1)), c("N", paste0("N_C", seq_along(x)))) })
      cell_counts = cell_counts %>% reshape::melt(measure.vars=c())
      wh_cols = c("L1","L2", grep("N", colnames(cell_counts), value = TRUE))
      cell_counts = cell_counts[, wh_cols, drop=FALSE]
      colnames(cell_counts)[1:2] = c("j","k")
      
      stopifnot(all(cell_counts$j == sim_seeds$j & cell_counts$k == sim_seeds$k))
      realisation_data = cbind(sim_seeds, cell_counts[,-c(1:2)])
      realisation_data = realisation_data[!is.na(realisation_data$seed),]
      
      
      # observation data
      sample_list = lapply(particle_obs_list, lapply, lapply, "[[", "samples")
      dist_list = lapply(data$particles, "[[", "distances")
      
      sample_seeds = lapply(sample_list, lapply, lapply, "[[", "seed")
      sample_seeds = lapply(sample_seeds, lapply, lapply, null_to_na)
      
      sample_mrate_sf = lapply(particle_obs_list, lapply, lapply, "[[", "mutation_rate")
      sample_mrate_sf = lapply(sample_mrate_sf, lapply, lapply, null_to_na)
      
      l_level1 = length(sample_seeds)
      l_level2 = unique(unlist(lapply(sample_seeds, length)))
      l_level3 = unique(unlist(lapply(sample_seeds, lapply, length)))
      stopifnot(length(l_level2) == 1 & length(l_level3) == 1)
      
      idx_j = rep(seq_len(l_level1), each=l_level2*l_level3)
      idx_k = rep(rep(seq_len(l_level2), each=l_level3), l_level1)
      idx_l = rep(rep(seq_len(l_level3), l_level2), l_level1)
      
      idx_list = split(cbind(idx_j, idx_k, idx_l), seq_along(idx_l))
      dists = sapply(idx_list, function(x) { dist_list[[x]] })
      
      observation_data = data.frame(
        j = idx_j,
        k = idx_k,
        l = idx_l,
        sample_seed = unlist(sample_seeds),
        mutation_rate_scale_factor = unlist(sample_mrate_sf),
        dist = dists
      )
      
      
      # final list
      res_list = list(
        states = tibble::as_tibble(state_data),
        particles = tibble::as_tibble(particle_data),
        realisations = tibble::as_tibble(realisation_data),
        observations = tibble::as_tibble(observation_data)
      )
      
      return(res_list)
    }
  
  collapse_chain_list = function(x) {
    # collapse smc chain to one large tibble
    chain_el = c("states","particles","realisations","observations")
    chain_lists = rep(list(), length(chain_el))
    
    # list of element column names 
    chain_el_cns = 
      lapply(chain_el, function(e) {
        lapply(seq_along(x), function(i) colnames(x[[i]][[e]])) %>% 
          unlist() %>% unique()
      }) %>% magrittr::set_names(chain_el)
    
    for (i in seq_along(x)) {
      for (el in chain_el) {
        
        el_cns = chain_el_cns[[el]]
        el_d_i = x[[i]][[el]]
        
        # add na values for missing columns
        missing_el_cns = el_cns[!el_cns %in% colnames(el_d_i)]
        col_data_add = rep(NA, length(missing_el_cns))
        names(col_data_add) = missing_el_cns
        el_d_i = tibble::add_column(el_d_i, !!!col_data_add)
        
        c_chain_el = tibble::add_column(el_d_i[,el_cns], i=i, .before=1)
        chain_lists[[el]] = rbind(chain_lists[[el]], c_chain_el)
      }
    }
    return(chain_lists)
  }
  
  
  # return nothing if less than two states available
  particle_files = ordered_particle_list(dir)
  if (length(particle_files) < 2) return(NULL)
  
  
  # variables and target
  variables = file.path(dir, "smc_chain", "options.rds") %>% readRDS()
  target_tree = variables$target_y
  final_epsilon = readRDS(tail(particle_files, 1))$epsilon
  
  
  # extract ground truth values:
  case_dir = basename(dirname(dir))
  encodes_gt_data = grepl("-[[:alpha:]]_[[:digit:].]+", dir)
  
  if (encodes_gt_data) {
    
    string_el = strsplit(strsplit(case_dir, "-")[[1]], "_")
    wh_param_el = sapply(string_el, length) == 2
    
    param_values = 
      sapply(string_el[wh_param_el], "[[", 2) %>% 
      magrittr::set_names(sapply(string_el[wh_param_el], "[[", 1))

    ids = c(
      N = NA,
      mutationrate = "mutation_rates.1",
      push = "push_power.1",
      death = "deathrate.1",
      selection = "birthrates.2",
      time = "clone_start_times.2"
    )[names(param_values)]
    
    data_gt =
      data.frame(value = as.numeric(param_values), parameter = ids) %>% 
      dplyr::filter(!is.na(parameter))

  } else {
    data_gt = NULL
  }
  

  # load smc chain
  cat("-> Loading SMC chain\n")
  chain = particle_files %>% 
    pbmcapply::pbmclapply(load_particle_set, eps=final_epsilon, mc.cores=n_cores) %>% 
    (function(x) x[!sapply(x, is.null)]) %>%
    collapse_chain_list()
  cat("\n")
  
  
  # add labels to chain data:
  for (el in names(chain)) {
    if ("i" %in% colnames(chain[[el]])) {
      id = names(particle_files)[unlist(chain[[el]][,"i"])]
      id = factor(id, names(particle_files), ordered = TRUE)
      chain[[el]] = tibble::add_column(chain[[el]], i_id=id, .before=1)
    }
  }
  
  # mark duplicated realizations
  wh_cols = colnames(chain$realisations)[!colnames(chain$realisations) %in% c("i_id","i","j","k")]
  chain$realisations$duplicated = duplicated(chain$realisations[,wh_cols])
  
  # mark duplicated observations
  wh_dup = chain$realisations$duplicated
  idx_dup = with(chain$realisations[wh_dup,], paste0(i, ",", j, ",", k))
  idx_obs =  with(chain$observations, paste0(i, ",", j, ",", k))
  chain$observations$duplicated = idx_obs %in% idx_dup
  
  
  # load rejected particle chain
  if (include_rejected_particles) {
    cat("-> Loading rejected particle chain\n")
    
    rejected_particle_files = ordered_particle_list_rejected(dir)
    
    # load smc chain
    cat("-> Loading SMC chain\n")
    chain_rejected = rejected_particle_files %>% 
      pbmcapply::pbmclapply(load_particle_set, eps=final_epsilon, mc.cores=n_cores) %>% 
      (function(x) x[!sapply(x, is.null)]) %>%
      collapse_chain_list()
    cat("\n")
    
    
    # add labels to chain data:
    for (el in names(chain_rejected)) {
      if ("i" %in% colnames(chain_rejected[[el]])) {
        id = names(particle_files)[unlist(chain_rejected[[el]][,"i"])]
        id = factor(id, names(particle_files), ordered = TRUE)
        chain_rejected[[el]] = tibble::add_column(chain_rejected[[el]], i_id=id, .before=1)
      }
    }
    
    names(chain_rejected) = paste0("rejected-", names(chain_rejected))
    
    cat("\n")
  } else {
    chain_rejected = NULL
  }
  
  result = c(list(variables=variables, target=target_tree, ground_truth=data_gt), chain, chain_rejected)
  attr(result, "class") = "chess_result_set"
  
  cat("-> Done.\n")
  return(result)
  
}


update_smc_result_set = function(sim_dir, force=FALSE, ...) {
  
  # indentify chain files
  dump_file = file.path(sim_dir, "result_set.rds")
  mtime_particle_files = ordered_particle_list(sim_dir) %>% file.mtime()
  mtime_last_particle_file = max(mtime_particle_files)
  
  
  # check age of dump file
  if (file.exists(dump_file)) { 
    mtime_dump_file = file.mtime(dump_file)
    new_states = any(mtime_dump_file < mtime_particle_files)
  } else {
    new_states = TRUE
  }
  
  
  if (!file.exists(dump_file) | force | new_states) { 
    
    results_this = 
      tryCatch({
        load_smc_result_set(sim_dir, ...)
      }, error=function(e) {
        print(e)
        return(NULL)
      })
    
    if (!is.null(results_this)) {
      saveRDS(results_this, dump_file, version=2)
      Sys.setFileTime(dump_file, mtime_last_particle_file)
    }
    
  } else { # (re)load the result set from saved chains
    results_this = dump_file
  }
  
  invisible(results_this)
}

ordered_particle_list = function(dir) {
  
  chain = dir %>% 
    list.files("smc_chain", full.names=1) %>%
    list.files("state_step", full.names=1)
  
  if (length(chain) == 0) return(chain)
  
  # set names 
  state_num = as.numeric(sub("[.]rds$", "", gsub(".*_", "", basename(chain))))
  chain_num = as.numeric(gsub("smc_chain", "0", basename(dirname(chain))))
  names(chain) = paste0("S", chain_num, ".", state_num)
  
  chain[order(chain_num, state_num)] # order
}


ordered_particle_list_rejected = function(dir) {
  
  chain = dir %>% 
    list.files("rejected", full.names = TRUE, recursive = TRUE, include.dirs = TRUE) %>%
    list.files("S", full.names=TRUE)
  
  if (length(chain) == 0) return(chain)
  
  # set names 
  state_num = as.numeric(gsub("[.]rds$", "", gsub("S", "", basename(chain))))
  chain_num = as.numeric(gsub("smc_chain", "0",basename(dirname(dirname(chain)))))
  names(chain) = paste0("S", chain_num, ".", state_num)
  
  chain[order(chain_num, state_num)] # order
  
}

plot_smc_abc_result = function(x, plot_chain=FALSE, title="") {
  
  stopifnot("chess_result_set" %in% class(x))
  
  data_long = 
    x$particles %>% 
    get_varied_parameter_list() %>% 
    dplyr::select(-frac_below_eps) %>%
    reshape::melt(id.vars=c("i","weight","weight_i"))
  
  if (!plot_chain) {
    
    u_states = sort(unique(data_long$i))
    
    data_long$i = factor(
      data_long$i,
      levels = c(min(u_states), max(u_states)),
      labels = c("prior", "posterior"),
      ordered = TRUE
    )
    
  } else {
    
    data_long$i = factor(
      data_long$i,
      x$states$i, 
      x$states$i_id,
      ordered = TRUE
    )
  }
  
  
  posterior_mode = 
    split(data_long, data_long$i) %>% 
    lapply(function(d) split(d, d$variable)) %>% 
    lapply(lapply, function(d) {wh = !is.na(d$value); if (sum(wh) == 0) return(NA); de=density(d$value[wh], weights=d$weight_i[wh]); return(de$x[which.max(de$y)])}) %>% 
    #lapply(lapply, function(d) weighted.mean(d$value, d$weight)) %>%
    reshape2::melt() %>%
    magrittr::set_colnames(c("value","variable","i"))
  
  
  plot_posterior = 
    data_long %>%
    dplyr::filter(!is.na(i)) %>% 
    #mutate(type=factor(type, c("prior","posterior"), ordered = TRUE)) %>%
    ggplot(aes(y=value, x=i)) + 
    facet_wrap(~variable, scales="free", ncol=1) + 
    geom_violin(aes(weight=weight_i), scale="width", fill="gray40", width=0.5) + 
    geom_boxplot(aes(weight=weight_i), width=0.1, fill="gray90", color="gray5", outlier.shape = NA) + 
    geom_point(data=posterior_mode, color="gray5", size=1) + 
    geom_line(data=posterior_mode, aes(group=1), color="gray5", linetype=2) +
    xlab("") + 
    ylab("Value") + 
    ggtitle(title) + 
    theme(axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5))
  
  
  return(plot_posterior)
}

plot_smc_abc_scatter = function(x, title="") {
  dwide = get_varied_parameter_list(x)
  dwide$weight = dwide$frac_below_eps / sum(dwide$frac_below_eps, na.rm = TRUE)
  dwide$weight_ = dwide$weight * runif(NROW(dwide), 0.975, 1.025)
  plot(data.frame(dwide[,!colnames(dwide) %in% c("i","weight","weight_i","frac_below_eps")]), cex=dwide$weight/max(dwide$weight) + 0.2)
  invisible(dwide)
}

plot_smc_abc_eps = function(data, title="") {
  
  state_data = data$states
  
  # add first observation distance if infinite:
  wh_first = state_data$i == 1
  if (is.infinite(state_data$epsilon[wh_first])) {
    state_data$epsilon[wh_first] = 
      max(data$observations$dist[data$observations$i == 1])
  }
  
  if (!"i_id" %in% names(state_data)) {
    state_data$i_id = factor(
      as.character(state_data$i),
      as.character(state_data$i),
      paste0("S", seq_along(state_data$i) - 1),
      ordered = TRUE
    )
  }
  
  state_data %>% 
    dplyr::mutate(epsilon=ifelse(i==1, NA, epsilon)) %>% 
    ggplot(aes(x=i_id, y=epsilon)) +
    geom_line(aes(group=1), linetype=2) +
    geom_point() +
    xlab("State") + 
    ylab(bquote(epsilon)) + 
    scale_y_continuous(labels = scales::scientific)
  
}

plot_smc_abc_eps_delta = function(data, title="") {
  
  state_data = data$states
  
  # add first observation distance if infinite:
  wh_first = state_data$i == 1
  if (is.infinite(state_data$epsilon[wh_first])) {
    state_data$epsilon[wh_first] = 
      max(data$observations$dist[data$observations$i == 1])
  }
  
  
  if (!"i_id" %in% names(state_data)) {
    state_data$i_id = factor(
      as.character(state_data$i),
      as.character(state_data$i),
      paste0("S", seq_along(state_data$i) - 1),
      ordered = TRUE
    )
  }
  
  # add delta
  d_eps = c(NA, diff(state_data$epsilon))
  eps_n_prev = c(NA, state_data$epsilon[-NROW(state_data)])
  state_data$delta_epsilon = -d_eps

  delta_epsilon_min = data$variables$delta_epsilon_min
  if (is.null(delta_epsilon_min)) { delta_epsilon_min = 0 }
  deleta_min = delta_epsilon_min  / 10
  state_data$rel_delta_epsilon = -d_eps / eps_n_prev + deleta_min
  state_data$below_delta = state_data$rel_delta_epsilon <= (delta_epsilon_min + deleta_min)
  state_data$below_delta[is.na(state_data$below_delta)] = FALSE
  
  state_data %>% 
    ggplot(aes(x=i_id, y=rel_delta_epsilon)) +
    geom_line(aes(group=1), linetype=2) +
    geom_point(aes(color=below_delta)) +
    xlab("State") + 
    ylab(bquote((epsilon[i-1]-epsilon[i])/epsilon[i-1])) + 
    geom_hline(yintercept=delta_epsilon_min, linetype=3) + 
    scale_color_manual(breaks=c("TRUE","FALSE"), values=c("gray30","black")) + 
    guides(color=FALSE)  + 
    scale_y_log10()
    
}


plot_smc_abc_ess = function(data, title="") {
  
  state_data = data$states

  if (!"i_id" %in% names(state_data)) {
    state_data$i_id = factor(
      as.character(state_data$i),
      as.character(state_data$i),
      paste0("S", seq_along(state_data$i) - 1),
      ordered = TRUE
    )
  }
  
  state_data %>%
    dplyr::mutate(ess=ifelse(i==1, NA, ess)) %>% 
    ggplot(aes(x=i_id, y=ess)) +
    geom_line(aes(group=1), linetype=2) +
    geom_point() +
    xlab("State") + 
    ylab("ESS")  +  
    geom_hline(yintercept=data$variables$N_T, linetype=3) + 
    ylim(0, data$variables$N)
  
}

plot_smc_abc_acceptance = function(data, title="") {
  
  state_data = data$states
  
  if (!"i_id" %in% names(state_data)) {
    state_data$i_id = factor(
      as.character(state_data$i),
      as.character(state_data$i),
      paste0("S", seq_along(state_data$i) - 1),
      ordered = TRUE
    )
  }
  
  state_data %>% 
    ggplot(aes(x=i_id, y=acceptance_rate)) +
    geom_line(aes(group=1), linetype=2) +
    geom_point() +
    xlab("State") + 
    #ylab("Acceptance rate") + 
    ylim(0, 1)  +
    geom_hline(yintercept=data$variables$min_acceptance_rate, linetype=3) +
    ylab("Accept. rate")
  
}

plot_smc_abc_states = function(res, title="", include_delta_eps=TRUE) {
  
  blank_x =
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  p_ess =  plot_smc_abc_ess(res) + blank_x
  p_eps = plot_smc_abc_eps(res) + blank_x 
  
  if (include_delta_eps) {
    p_eps_delta = plot_smc_abc_eps_delta(res) +  blank_x
  } else{
    p_eps_delta = NULL
  } 
    
  p_acc = plot_smc_abc_acceptance(res) 
  p_acc = p_acc + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p_acc = p_acc + labs(caption = paste0("Final tolerance: ", signif(tail(res$states$epsilon, 1), 5)))
  
  
  plot_states = 
    (gtable_rbind(
      ggplotGrob(p_ess),
      ggplotGrob(p_eps),
      ggplotGrob(p_eps_delta),
      ggplotGrob(p_acc),
      size = "max"
    ))
  
}

plot_abc_result = function(particle_list, prior_dist=NULL, title="") {
  
  short_to_long = function(d) {
    d %>% 
      lapply(function(x) x$particle$param.matrix) %>% 
      reshape::melt() %>% 
      dplyr::mutate(parameter=paste0(X1, ".", X2)) %>% 
      select(-X1, -X2)
  }
  
  
  params = reshape::melt(lapply(particle_list, short_to_long))
  varied = tapply(params$value, params$parameter, function(x) any(x!=x[1]))
  varied_ids = names(which(varied))
  
  if (!is.null(prior_dist)) {
    prior_dist = 
      short_to_long(prior_dist)  %>%  
      dplyr::filter(parameter %in% variable_id) %>%
      dplyr::mutate(type="prior") %>%
      dplyr::select(parameter, type, value)
  }
  
  post_dist = 
    particle_list %>% 
    dplyr::filter(parameter %in% variable_id) %>% 
    dplyr::mutate(type="posterior") %>% 
    dplyr::select(parameter, type, value)
  
  posterior_mean = 
    tapply(post_dist$value, post_dist$parameter, function(d) {de=density(d); return(de$x[which.max(de$y)])}) %>% 
    reshape2::melt() %>%
    magrittr::set_colnames(c("parameter","value")) %>% 
    dplyr::mutate(type="posterior")
  
  plot_posterior = 
    rbind(prior_dist, post_dist) %>%
    dplyr::mutate(type=factor(type, c("prior","posterior"), ordered = TRUE)) %>%
    ggplot(aes(y=value, x=type)) + 
    facet_wrap(~parameter, scales="free") + 
    geom_violin(scale="width", fill="gray40", width=0.5) + 
    geom_boxplot(width=0.1, fill="gray90", color="gray5") + 
    geom_point(data=posterior_mean, color="gray5", size=1) + 
    xlab("Parameter") + 
    ylab("Value") + 
    ggtitle(title)
  
  return(plot_posterior)
}

head.chess_result_set = function(x, n=10, ...) {
  x$observations = head(dplyr::arrange(x$observations, dist), n)
  drop_unobserved(x)
}

tail.chess_result_set = function(x, n=10, ...) {
  x$observations = tail(dplyr::arrange(x$observations, dist), n)
  drop_unobserved(x)
}

drop_unobserved = function(x) {
  
  idx_obs = apply(x$observations[,c("i","j","k")], 1, paste0, collapse = ",")
  idx_real = apply(x$realisations[,c("i","j","k")], 1, paste0, collapse = ",")
  
  x$realisations = x$realisations[idx_real %in% idx_obs, ]
  idx_real = apply(x$realisations[,c("i","j")], 1, paste0, collapse = ",")
  idx_part = apply(x$particles[,c("i","j")], 1, paste0, collapse = ",")
  
  x$particles = x$particles[idx_part %in% idx_real,]
  
  return(x)
}

plot_smc_results = function(input, dir=NULL, n_top_sims=25, n_top_trees=50, map_distances=FALSE, collect_tree_set=0, n_cores=1, labeller_function=function(x) return(x), plot_tree_function=CHESS::plot_tree, rescale_generation_tree=TRUE, n_clonal=NULL, run_mobster=FALSE, add_history_plot=TRUE, dump_data=FALSE, ...) {
  
  # either load chain from disk if input is a directory 
  # or assume it is a pre-loaded chain
  if (is.character(input)) {
    if (file.exists(input)) {
      result_set = readRDS(input)
    } else {
      result_set = load_smc_result_set(dir)
    }
  } else {
    stopifnot(!is.null(dir))
    result_set = input
  }
  
  # overload ggsave function if plot objects should be dumped
  if (dump_data) {
    ggsave = function(filename, plot, ...) {
      of = file.path(dirname(filename), paste0(".", gsub("[.][a-zA-Z]+$", ".rds", basename(filename))))
      saveRDS(plot, of)
      ggplot2::ggsave(filename, plot, ...)
    }
  }
  
  stopifnot("chess_result_set" %in% class(result_set))
  
  # temporary fix:
  if (!is.null(result_set$ground_truth)) {
    wh = result_set$ground_truth$parameter == "mutation_rates.1"
    tmp = result_set$ground_truth[wh,]
    tmp$parameter = "mutation_rate"
    result_set$ground_truth = rbind(result_set$ground_truth, tmp)
    
    result_set$ground_truth = 
      result_set$ground_truth %>% 
      dplyr::mutate(parameter = gsub("deathrate[.]", "deathrates.", parameter))
  }
  
  # add frac below epsilion to data:
  eps = min(result_set$states$epsilon)
  result_set = recalculate_fracs_below_eps(result_set, eps)
  
  #####################################
  
  if (map_distances) {
    
    weight = sapply(result_set$chain$S1$particles, function(x) x$density)
    params = lapply(result_set$chain$S1$particles, function(x) x$param.matrix)
    params_mean = apply(abind::abind(params, along = 3), c(1,2), weighted.mean, w=weight)
    const_params = c(result_set$variables$const_params, list(param.matrix=params_mean))
    seeds = sample(1e9, 10, replace = TRUE)
    
    particles = lapply(seeds, function(s) c(const_params, list(seed=s)))
    sgen = result_set$variables$s_gen
    msim = result_set$variables$M_sim
    obs_particles = observe_particles(particles, sgen, M=1, M_sim = msim, sample_tree = result_set$variables$sample_tree)
    
    rho = result_set$variables$rho
    mean_rho = list()
    for (i in seq_along(obs_particles)) {
      print(i)
      mean_rho[[i]] = numeric()
      #for (j in seq_along(obs_particles[[i]]$observations)) {
      j=1
      for (k in seq_along(obs_particles[[i]]$observations[[j]])[1:10]) {
        target_y = obs_particles[[i]]$observations[[j]][[k]]$tree
        d_matrix = d_matrix_from_particles(obs_particles[-i], target_y, rho)
        mean_rho[[i]] = c(mean_rho[[i]], c(d_matrix))
      }
      #}
    }
    
    plot_dist_distribution = 
      mean_rho %>% 
      reshape2::melt() %>% 
      ggplot(aes(x=factor(L1, 1:100), y=value)) + 
      geom_violin(scale="width", fill="gray40", width=0.75) + 
      geom_boxplot(width=0.2, fill="gray90", color="gray5", outlier.shape = NA) +
      #geom_point(data=posterior_mean, aes(y=value, x=L1), inherit.aes = FALSE, color="gray5", size=1) + 
      #geom_line(data=posterior_mean, aes(y=value, x=L1, group=1), inherit.aes = FALSE, color="gray5", linetype=2) +
      xlab("Sim #") + 
      ylab("Expected distance") + 
      geom_hline( aes(linetype=factor("Achieved distance to target"), yintercept = result_set$chain$S1$epsilon)) + 
      scale_linetype_manual(values = 2, breaks="Achieved distance to target") + 
      labs(color="", linetype="") +
      theme(legend.position="bottom") 
    
    out_file = file.path(dir, "posterior_mean_distance_distribution.pdf")
    ggsave(out_file, plot_dist_distribution, height=2.5, width=4)
    
    dist_data_save = list(distances=mean_rho, archived=result_set$chain$S1$epsilon)
    out_file = file.path(dir, "posterior_mean_distance_distribution.rds")
    saveRDS(dist_data_save, out_file, version = 2)
    
  }
  
  ####################################
  
  try({
    
    plot_variables = plot_smc_abc_result(result_set, plot_chain=T) 
    
    if (!is.null(result_set$ground_truth)) {
      plot_variables =
        plot_variables +
        geom_hline(
          data = result_set$ground_truth %>%
            dplyr::filter(parameter %in% plot_variables$data$variable) %>%
            dplyr::mutate(variable = parameter),
          aes(yintercept = value),
          linetype = 3
        )
    }
    
    height = length(unique(plot_variables$data$variable)) * 1.5 + 1.5
    out_file = file.path(dir, "smc_abc_variable_histograms.pdf")
    ggsave(out_file, plot_variables, height=height)
  })
  
  ####################################
  
  try({
    
    plot_variables_post = plot_smc_abc_result(result_set)
    
    if (!is.null(result_set$ground_truth)) {
      plot_variables_post =
        plot_variables_post +
        geom_hline(
          data = result_set$ground_truth %>%
            dplyr::filter(parameter %in% plot_variables_post$data$variable) %>%
            dplyr::mutate(variable = parameter),
          aes(yintercept = value),
          linetype = 3
        )
    }
    
    height = length(unique(plot_variables_post$data$variable)) * 1.5 + 1.5
    out_file = file.path(dir, "smc_abc_variable_histograms_posterior.pdf")
    ggsave(out_file, plot_variables_post, width = 2.5, height = height)
    
  })
  
  ####################################
  
  try({
    plot_states = plot_smc_abc_states(result_set)
    out_file = file.path(dir, "smc_abc_states.pdf")
    n_states = NROW(result_set$states)
    ggsave(out_file, plot_states, width=0.95+(3.05*n_states/18), height=5.8)
  })
  
  #####################################
  
  tryCatch({
    pdf(file.path(dir, "pairwise_scatter.pdf"), width = 7, height=7)
    result_set$particles %>% dplyr::filter(i==max(i)) %>% plot_smc_abc_scatter()
  }, finally = {
    dev.off()
  }, error = function(e) print(e))
  
  tryCatch({
    pdf(file.path(dir, "pairwise_scatter_all_states.pdf"), width = 7, height=7)
    result_set$particles %>% plot_smc_abc_scatter()
  }, finally = {
    dev.off()
  }, error = function(e) print(e))
  
  #####################################
  
  
  try({
    
    ppd =
      result_set$particles %>%
      dplyr::filter(i == max(i)) %>%
      posterior_density_plot(result_set$particles, gt = result_set$ground_truth)
    
    plot_width = 3.2 * length(ppd$widths) + 0.5
    plot_height = 3.0 * length(ppd$heights)
    out_file = file.path(dir, "posterior_density.png")
    ggsave(out_file, ppd, width = plot_width, height = plot_height)
    
  })  
  
  try({
    ppd = posterior_density_plot(result_set$particles, gt = result_set$ground_truth)
    plot_width = 3.2 * length(ppd$widths) + 0.5
    plot_height = 3.0 * length(ppd$heights)
    out_file = file.path(dir, "posterior_density_all_states.png")
    ggsave(out_file, ppd, width = plot_width, height = plot_height)
  })
  
  ###################################

  target = result_set$variables$target_y
  rho = result_set$variables$rho
  tree_distance_mod = function(...) tree_distance(..., return_tree=TRUE)
  assign("tree_distance", tree_distance_mod, envir=environment(rho)) # patch distance function
  
  modify_tree_function =
    function(x, ...) {
      if ("tree" %in% names(x)) {
        x$tree =  rho(target, x$tree)$tree
      } else {
        x = rho(target, x)$tree
      }
      return(x)
    }
  
  ###################################
  
  sorted_obs_list = 
    dplyr::arrange(result_set$observations, dist) %>% 
    dplyr::filter(!duplicated)
  
  idx = unlist(data.frame(sorted_obs_list[1,][c("i","j","k","l")]))
  
  try({
    
    plot_idx =
      replot_particle_by_idx(
        result_set,
        idx,
        alpha_burden = FALSE,
        target_tree = result_set$variables$target_y,
        labeller_function = labeller_function,
        sample_tree = result_set$variables$sample_tree,
        plot_tree_function = plot_tree_function,
        modify_tree_function = modify_tree_function,
        rescale_generation_tree = rescale_generation_tree,
        n_clonal = n_clonal,
        run_mobster = run_mobster,
        add_history_plot = add_history_plot,
        ...
      )
  
    out_file = file.path(dir, "best_fit.pdf")
    width = ifelse(run_mobster, 15, 12.65) - ifelse(add_history_plot, 0, 3.4)
    ggsave(out_file , plot_idx, width = width, height=3.3, bg = "transparent")
    
  })
  
  ####################################

  if (n_top_trees) {
    try({
      n_top_trees = min(c(n_top_trees, NROW(sorted_obs_list)))
      
      # plot 50 closes trees into one pdf
      cat(paste0("Plotting ", n_top_trees, " closest trees.\n"))
      plot_list =
        pbmcapply::pbmclapply(seq_len(n_top_trees), function(i) {
          idx = unlist(data.frame(sorted_obs_list[i, ][c("i", "j", "k", "l")]))
          sim = get_sim_by_idx(result_set, idx, TRUE)
          tm = result_set$variables$tree_manipulator
          tree = get_tree(sim, x$variables$sample_tree, tree_manipulator=tm)
          tree_dist = rho(target, tree)
          if (!isTRUE(all.equal(sorted_obs_list$dist[i], as.numeric(tree_dist$dist))))
            stop("IMPORTANT: Can't recreate tree.\n")
          tree = remove_root_tip(tree_dist$tree, "GL")
          title = paste0(i, "/", n_top_trees)
          subtitle = paste0("d = ", signif(sorted_obs_list$dist[i], 3))
          tree_plot = plot_tree_function(tree) + ggtitle(title, subtitle)
          tree_plot = add_driver_labels_to_tree_plot(tree_plot, tree, sim)
          ggplotify::as.grob(tree_plot)
        }, mc.cores = n_cores, mc.preschedule = FALSE)
      
      grid_plot = gridExtra::arrangeGrob(grobs = plot_list, nrow = 10, ncol = 5)
      out_file = file.path(dir, "top_trees.pdf")
      ggsave(out_file, grid_plot, width = 10, height = 22)
    })
  }

  ####################################
  
  # plot the 25 closest simulations in full
  if (n_top_sims) {
    try({
      n_top_sims = min(c(n_top_sims, NROW(sorted_obs_list)))
      
      cat(paste0("Plotting ", n_top_sims, " closest simulations.\n"))
      
      plot_closest_n =
        pbmcapply::pbmclapply(seq_len(n_top_sims), function(i) {
          replot_particle_by_idx(
            result_set,
            unlist(data.frame(sorted_obs_list[i, ][c("i", "j", "k", "l")])),
            ,
            alpha_burden = FALSE,
            target_tree = result_set$variables$target_y,
            labeller_function = labeller_function,
            sample_tree = result_set$variables$sample_tree,
            plot_tree_function = plot_tree_function,
            modify_tree_function = modify_tree_function,
            rescale_generation_tree = rescale_generation_tree,
            n_clonal = n_clonal,
            run_mobster = run_mobster,
            add_history_plot = add_history_plot,
            ...
          ) %>% ggplotify::as.grob()
          
        }, mc.cores = n_cores, mc.preschedule = FALSE)
      
      grid_plot = gridExtra::marrangeGrob(grobs = plot_closest_n,
                                          nrow = 1,
                                          ncol = 1)
      out_file = file.path(dir, "top_sims.pdf")
      width = ifelse(run_mobster, 15, 12.65) - ifelse(add_history_plot, 0, 3.4)
      ggsave(out_file, grid_plot,  width = width, height = 3.5)
    })
  }
  
  #####################################
  
  if (collect_tree_set) {
    try({
      cat(paste0("Collecting ", collect_tree_set, " closest trees.\n"))
      
      trees = 
        pbmcapply::pbmclapply(seq_len(collect_tree_set), function(i) {
          get_tree_by_idx(
            result_set, 
            unlist(data.frame(sorted_obs_list[i,][c("i","j","k","l")]))
          ) %>% modify_tree_function()
        }, mc.cores = n_cores, mc.preschedule = FALSE)
       
      class(trees) = "multiPhylo"
      saveRDS(trees, file.path(dir, "posterior_trees.rds"))
      
    })
  }
  
  #####################################
}

dens_plot = function(data, a, b, lims=NULL) {
  
  wh = apply(is.na(data[,c(a,b)]), 1, any)
  x = data[!wh,a]
  y = data[!wh,b]
  w = data$weight[!wh]

  if (is.null(lims)) lims = c(range(x, na.rm=TRUE), range(y, na.rm=TRUE))
  kde = ggtern::kde2d.weighted(x, y, n = 1000, w = w, lims = lims)
  data_density = reshape::melt(kde$z) %>% dplyr::mutate(X1 = kde$x[X1], X2 = kde$y[X2])
  colnames(data_density) = c(a, b, "d")
  data_density$d = data_density$d / sum(data_density$d)
  
  plot =
    ggplot(data_density, aes(x=get(a), y=get(b))) +
    geom_tile(aes(fill=d)) +
    geom_contour(aes(z=d), alpha=0.75, color="gray5", linetype=2, bins=5) +
    #scale_fill_viridis_c(direction = 1) +
    scale_fill_distiller(palette = 9, direction = 1) +
    geom_point(data=data_density %>% dplyr::filter(d == max(d)), size=1, color="red")
  
  return(plot)
}

posterior_density_plot = function(x, x_original=NULL, gt=NULL) {
  
  data_wide = get_varied_parameter_list(x) 
  data_wide$weight = data_wide$frac_below_eps / sum(data_wide$frac_below_eps, na.rm = TRUE)
  
  if (!is.null(x_original)) {
    data_wide_prior = get_varied_parameter_list(x_original)
  } else {
    data_wide_prior = data_wide
  }
  
  data_wide = data_wide %>% dplyr::filter(weight > 0)
  
  min_values_prior = apply(data_wide_prior, 2, min, na.rm=TRUE)
  max_values_prior = apply(data_wide_prior, 2, max, na.rm=TRUE)

  cn = colnames(data_wide)[!colnames(data_wide) %in% c("weight","weight_i","i","frac_below_eps")]
  plot_list = list()
  
  for (cn2 in cn) {
    for (cn1 in cn) {
      
      if (cn1 == cn2) {
        
        plot = 
          ggplot(data_wide, aes(x=get(cn1))) + 
          geom_histogram(aes(y = ..density.., weight = weight), bins=100) + 
          xlab(cn1)
        
        if (!is.null(data_wide_prior)) {
          plot = plot + 
            xlim(min_values_prior[cn1], max_values_prior[cn1])
        }
        
        if (cn1 %in% gt$parameter) {
          plot = plot + 
            geom_vline(data=gt %>% filter(parameter == cn1), 
                       aes(xintercept=value), color="gray10", linetype=2)
        }
        
      } else if (match(cn1, cn) > match(cn2, cn)) {
        
        plot = 
          ggplot(data_wide, aes(x=get(cn1), y=get(cn2), size=weight)) +
          geom_point(alpha=0.25) + 
          xlab(cn1) + 
          ylab(cn2) + 
          guides(color=FALSE, size=FALSE) + 
          scale_color_viridis_c() + 
          scale_size(range =c(0.1, 1.5))
        
        if (!is.null(data_wide_prior)) {
          plot = plot + 
            xlim(0, max_values_prior[cn1]) + 
            ylim(0, max_values_prior[cn2])
        }
        
        
      } else {
        
        if (!is.null(data_wide_prior)) {
          lims = c(
            0, 
            max_values_prior[cn1],
            0, 
            max_values_prior[cn2]
          )
        } else {
          lims = NULL
        }
        
        plot = 
          dens_plot(data_wide, cn1, cn2, lims) +
          xlab(cn1) + 
          ylab(cn2) + 
          guides(fill=FALSE)
        
        
        if (!is.null(gt)) {
          if (all(c(cn1,cn2) %in% gt$parameter)) {
            gt_data = reshape2::dcast(gt, .~parameter)
            plot = plot + geom_point(data=gt_data, fill=NA, color="gray5", shape=4, size=3)
          } 
        }
        
      }
      
      idx = length(plot_list) + 1
      plot_list[[idx]] = ggplotify::as.grob(plot)
    }
  }
  
  
  arrangeGrob(grobs=plot_list, ncol=length(cn), nrow=length(cn))
}

get_chain_length = function(x) {
  tryCatch({
    sum(basename(ordered_particle_list(x)) != "state_step_0.rds")
  }, error=function(e) return(NA))
}

smc_abc_finished = function(x) {
  tryCatch({
    file.exists(file.path(x, "done.txt"))
  }, error=function(e) return(NA))
}


get_varied_parameter_list = function(p_data, include_rejected=FALSE) {
  
  # catch call with full result set as argument
  if ("chess_result_set" %in% class(p_data)) {
    if (include_rejected & "rejected-particles" %in% names(p_data)) {
      p_data = 
        rbind(
          p_data$particles %>% dplyr::mutate(type="accepted"), 
          p_data$`rejected-particles`  %>% dplyr::mutate(type="rejected")
        )
    } else {
      p_data = p_data$particles %>% dplyr::mutate(type="accepted")
    }
  }
  
  # clones parameters
  not_clone_params = c("i_id","i","j","x","y","z","weight","verbose","type","frac_below_eps")
  clone_params = colnames(p_data)[!colnames(p_data) %in% not_clone_params]
  
  # find varied parameters
  is_varied = sapply(clone_params, function(x) NROW(unique(p_data[x])) > 1)
  varied_params = names(is_varied)[is_varied]
  
  # drop some:
  varied_params = varied_params[!varied_params %in% c("mutation_rate_sf")]
  
  # return data 
  if (include_rejected) {
    params_ret = unique(c(varied_params, "i", "weight", "weight_i", "type","frac_below_eps"))
  } else {
    params_ret = unique(c(varied_params, "i", "weight", "weight_i","frac_below_eps"))
  }
  
  data.frame(p_data[params_ret])
}



print.chess_result_set = function(x, ...) {
  
  head(x$observations)
  
  n_dim = apply(x$observations[,c("i","j","k","l")], 2, function(x) length(unique(x)))
  dim_msg = paste0(n_dim, collapse="x")
  sum_particles = NROW(x$observations)
  msg = paste0("SMC result set containing ", dim_msg, " (N=", sum_particles, ") trees\n")
  cat(msg)
  
  eps_final = signif(tail(x$states$epsilon, 1), 3)
  n_below = sum(x$observations$dist <= eps_final)
  frac_below = signif(n_below / NROW(x$observations) * 100, 3)
  msg = paste0(" => below final eps=", eps_final, " N=", n_below, " (", frac_below,"%) of trees\n")
  cat(msg)
  
  
}



random_posterior_tree_indices = function(x, n=10, seed=NULL) {
  
  stopifnot("chess_result_set" %in% class(x))
  stopifnot((is.numeric(seed) | is.null(seed)) & length(seed) == 1 )
  
  if (!is.null(seed)) set.seed(seed)
  
  # find trees that are part of the posterior set and not a duplicate
  min_eps = min(x$states$epsilon)
  is_post = x$observations$dist < min_eps
  is_dup = duplicated(x$observations[,c("j","k","l","sample_seed","dist")])
  
  # sample subset 
  wh_can_sample = which(is_post & !is_dup)
  n_get = min(c(n, length(wh_can_sample)))
  idx_sampled = sample(wh_can_sample, n_get, FALSE)
  
  # convert indixes to list 
  idxs = as.list(data.frame(t(x$observations[idx_sampled, c("i","j","k","l")])))
  names(idxs) = NULL
  
  return(idxs)
}


top_posterior_tree_indices = function(x, n=10) {
  
  stopifnot("chess_result_set" %in% class(x))
  
  # find trees that are part of the posterior set and not a duplicate
  min_eps = min(x$states$epsilon)
  obs = head(arrange(x$observations, dist), n*20)
  
  is_post = obs$dist < min_eps
  is_dup = duplicated(obs[,c("j","k","l","sample_seed","dist")])
  obs = obs[is_post & !is_dup,]
  
  obs_sampled = head(obs, n)
  
  # convert indixes to list 
  idxs = as.list(data.frame(t(obs_sampled[, c("i","j","k","l")])))
  names(idxs) = NULL
  
  return(idxs)
}



get_sim = function(x, idx, ...) {
  UseMethod("get_sim", x)
}

get_sample_setup = function(x, idx, ...) {
  UseMethod("get_sample_setup", x)
}

get_obs = function(x, idx, ...) {
  UseMethod("get_obs", x)
}

get_tree.chess_result_set = function(x, idx, ...) { 
  get_tree_by_idx(x, idx, ...)
}

get_sim.chess_result_set = function(x, idx, ...) { 
  get_sim_by_idx(x, idx, ...)
}

get_sample_setup.chess_result_set = function(x, idx, ...) { 
  get_sample_setup_by_idx(x, idx, ...)
}

get_param_matrix = function(x) {
  wh_cols = colnames(x)[grepl("[.]", colnames(x))]
  row_id = gsub("[.].*", "", wh_cols)
  col_id = gsub(".*[.]", "", wh_cols)
  
  split(x[,wh_cols], seq_len(NROW(x))) %>% 
    lapply(unlist) %>% 
    lapply(tapply, list(row_id, col_id), c)
} 


get_sim_params = function(res, idx) {
  wh_row = res$particles$i == idx[1] & res$particles$j == idx[2]
  c_param = res$particles[wh_row,]
  
  coords = as.list(unlist(c_param[, c("x","y","z")]))
  param_names = rownames(res$variables$p_gen()$param.matrix)
  params = get_param_matrix(c_param) %>% lapply(function(x) x[param_names,,drop=FALSE])
  seed = with(res$realisations, seed[i==idx[1] & j==idx[2] & k==idx[3]])
  
  
  # get other params
  wh_other = !grepl("[.]", colnames(c_param)) & 
             !grepl("weight", colnames(c_param)) & 
             !colnames(c_param) %in% c("x","y","z","mutation_rate_sf","i","j","i_id","mutation_rate","seed")

  other_params = as.list(c_param[, wh_other])
  
  c(coords, list(param.matrix=params[[1]], seed=seed), other_params)
}

get_sample_seed = function(x, idx) {
  with(x$observations, sample_seed[i==idx[1]&j==idx[2]&k==idx[3]&l==idx[4]])
}

get_obs_params = function(x, idx) {
  wh = with(x$observations, i==idx[1]&j==idx[2]&k==idx[3]&l==idx[4])
  x$observations[wh,]
}

replot_particle_by_idx = function(x, idx, rescale_generation_tree=TRUE, sample_tree=FALSE, ...) {
  
 sim = 
   tryCatch({
    # get simulation
    sim = get_sim_by_idx(x, idx, TRUE, rescale_generation_tree)
    }, error=function (e) { # maybe created with different code?
      x$variables$const_params = list(x$variables$const_params, list(alt_edge_finding = TRUE))
      sim = get_sim_by_idx(x, idx, TRUE)
    })
 
  x$variables$sample_tree = sample_tree
  tree = get_tree_by_idx(x, idx)
 
  # check correctness of trees
  # tm = x$variables$tree_manipulator
  # tree = get_tree(sim, x$variables$sample_tree, tree_manipulator=tm)
  # d_exp = with(x$observations, dist[i==idx[1]&j==idx[2]&k==idx[3]&l==idx[4]])
  # d_obs = x$variables$rho(tree$tree, x$variables$target_y)
  # if ("dist" %in% names(d_obs)) d_obs = d_obs$dist
  # stopifnot(d_obs == d_exp)
  
  plot_sim = plot.CHESS_S3(sim, sample_tree = tree$tree, ...)
  return(plot_sim)
}


get_sim_by_idx = function(x, idx, add_samples=FALSE, scale_tree=FALSE) {
  
  # get sim
  params = get_sim_params(x, idx)
  sim = do.call(new_simulation, params) 
  
  if ((sim$params$return_generation_tree & scale_tree) | add_samples) {
    # get and add sample setup
    sample_seed = get_sample_seed(x, idx)
    samples = x$variables$s_gen(sim, sample_seed)
    sim$samples = samples$setup
    sim$obs_params = as.list(get_obs_params(x, idx)[c("mutation_rate_scale_factor","sample_seed")])
  }
  
  if (sim$params$return_generation_tree & scale_tree) {
    #scale generation tree so data are equivalent to tree data
    
    # get node to mutation number tree:
    tree = 
      get_tree(
        x = sim, 
        x$variables$sample_tree, 
        tree_manipulator=x$variables$tree_manipulator, 
        include_node_labels=TRUE
      )
    
    if ("tree" %in% names(tree)) tree = tree$tree
    node_to_muts = tree$edge.length
    names(node_to_muts) = c(tree$tip.label, tree$node.label)[tree$edge[,2]]

    # scale the tree
    sf = sim$obs_params$mutation_rate_scale_factor
    if ("dispersion" %in% rownames(sim$params$param.matrix)) {
      disp =  sim$params$param.matrix["dispersion", 1]
    } else {
      disp = 0  
    }
    sim$sim$ScaleMutationTree(sf, node_to_muts, disp)

    # modify parameters to make changes clear:
    sim$params$return_generation_tree = FALSE
    sim$obs_params$mutation_rate_scale_factor = 1
  }
  
  if (add_samples) {
    # add samples
    sim$phi = samples$vars$phi
    sim = add_all_samples(sim, samples$setup)
  }
  
  return(sim)
  
}


get_sample_setup_by_idx = function(x, idx) {
  sim = get_sim_by_idx(x, idx)
  sample_seed = get_sample_seed(x, idx)
  x$variables$s_gen(sim, sample_seed)$setup
}

get_tree_by_idx = function(x, idx, use_alt_version=FALSE) {

  if (use_alt_version) {
    x$variables$const_params = list(x$variables$const_params, list(alt_edge_finding=TRUE))
  }
  
  sim = get_sim_by_idx(x, idx, TRUE)
  tm = x$variables$tree_manipulator
  
  tree = get_tree(sim, x$variables$sample_tree, tree_manipulator=tm)
  d = x$variables$rho(x$variables$target_y, tree$tree)
  if ("dist" %in% names(d)) d = d$dist
  dist = with(x$observations, dist[i==idx[1]&j==idx[2]&k==idx[3]&l==idx[4]])

  if (dist != d) { # try old/different version of lp assignment function.
    tree = get_tree(sim, x$variables$sample_tree, tree_manipulator=tm, insert_lp_above=TRUE)
    d = x$variables$rho(x$variables$target_y, tree$tree)
    if ("dist" %in% names(d)) d = d$dist
    stopifnot(dist == d)
  }
  
  return(tree)
} 

recalculate_fracs_below_eps = function(x, eps) {
  
  # recalc weight of particles
  idx_obs = with(x$observations, paste0(i, "-", j))
  frac = with(x$observations, tapply(dist <= eps, idx_obs, mean, na.rm=TRUE))
  
  # add weight
  idx_particles = with(x$particles, paste0(i, "-", j))
  x$particles$frac_below_eps = frac[idx_particles]
  
  if ("rejected-observations" %in% names(x)) {
    # recalc weight of particles
    idx_obs = with(x$`rejected-observations`, paste0(i, "-", j))
    frac = with(x$`rejected-observations`, tapply(dist <= eps, idx_obs, mean, na.rm=TRUE))
    
    # add weight
    idx_particles = with(x$`rejected-particles`, paste0(i, "-", j))
    x$`rejected-particles`$frac_below_eps = frac[idx_particles]
  }
  
  return(x)
}
