load_model_selection_data = function(dir, fname="model_selection_data_v3[.]rds") {
  files = list.files(dir, fname, recursive = FALSE, full.names = TRUE)
  files = files[!grepl("bad", files)]
  d = lapply(files, function(x) tryCatch(readRDS(x), error=function(e) return(NULL)))
  names(d) = basename(dirname(files))
  d = d[!sapply(d, is.null)]
  return(d)
}

plot_dic_values = function(d) {
  
  best_eps =
    d %>% 
    dplyr::mutate(val=pD-k) %>% 
    split(.$L1) %>% 
    lapply(function(x) x[which.min(abs(x$val)),]) %>% 
    do.call(what=rbind)
  
  best_eps2 =
    d %>% 
    dplyr::mutate(val=pD2-k) %>% 
    split(.$L1) %>% 
    lapply(function(x) x[which.min(abs(x$val)),]) %>% 
    do.call(what=rbind)
  
  avg_best_eps = mean(best_eps$eps, na.rm=TRUE)
  
  dic_closest =
    d %>% 
    split(.$L1) %>% 
    lapply(function(x) x[which.min(abs(x$eps - avg_best_eps)),]) %>% 
    do.call(what=rbind)
  
  plot_stat_values_across_eps(d, c("DIC","pD","DIC2","pD2")) + 
    geom_point(data=best_eps %>% dplyr::mutate(variable = "pD"), aes(y=pD)) + 
    geom_point(data=best_eps %>% dplyr::mutate(variable = "pD2"), aes(y=pD2)) + 
    geom_hline(data=best_eps %>% dplyr::mutate(variable = "pD"), aes(yintercept=k, color=L1), linetype=1, alpha=0.25) +
    geom_hline(data=best_eps2 %>% dplyr::mutate(variable = "pD2"), aes(yintercept=k, color=L1), linetype=1, alpha=0.25) +
    geom_vline(data=best_eps %>% dplyr::mutate(variable = "DIC"), aes(xintercept=eps, color=L1), linetype=1, alpha=0.25) + 
    geom_vline(data=best_eps2 %>% dplyr::mutate(variable = "DIC2"), aes(xintercept=eps, color=L1), linetype=1, alpha=0.25)
}

plot_stat_values_across_eps = function(d, stat) {
  
  d = d %>% 
    reshape2::melt(measure.vars=stat) %>% 
    dplyr::mutate(value=ifelse(L1 %in% c("DIC","DIC2"), log(value), value)) %>%
    dplyr::mutate(L1=ifelse(L1 %in% c("DIC","DIC2"), paste0("log(", L1, ")"), L1)) %>%
    dplyr::group_by(eps, L1, variable) %>% 
    dplyr::summarise(
      mean = mean(value), 
      low = quantile(value, 0.025), 
      high = quantile(value, 0.975),
      .groups = "drop"
    ) 
  
  d %>%
    arrange(eps) %>% 
    ggplot(aes(x=eps)) + 
    geom_ribbon(aes(ymin=low, ymax=high, fill=L1), alpha=0.3) + 
    geom_line(aes(y=mean, color=L1)) + 
    facet_wrap(~variable, scales="free") + 
    scale_color_brewer(palette = "Set1") + 
    scale_fill_brewer(palette = "Set1") + 
    labs(color="Model", fill="Model") + 
    xlab(bquote("Critical distance"~epsilon)) + 
    ylab("Value")
}

calc_dic_across_eps = function(d, n_cores=1, n_points=100, use_sl=FALSE) {
  calc_stat_across_eps(d, calc_dic, n_cores=n_cores, n_points=n_points, use_sl=use_sl, ...)
}

calc_aic_across_eps = function(d, n_cores=1, n_points=100, use_sl=FALSE, ...) {
  calc_stat_across_eps(d, calc_aic, n_cores=n_cores, n_points=n_points, use_sl=use_sl, ...)
}

calc_stat_across_eps = function(d, f, n_cores, n_points, ...) {
  
  min_eps = 0
  
  max_eps = 
    max(sapply(d, function(x)
      quantile(unlist(c(
        unlist(x$d_dhat), unlist(x$d_dhat_permutated)
      ))[!is.infinite(unlist(c(
        unlist(x$d_dhat), unlist(x$d_dhat_permutated)
      )))], 0.99, na.rm = TRUE))
    )
  
  eps_vals = seq(from=min_eps, to=max_eps, length=n_points)
  
  stats = 
    lapply(d, function(y) {
      do.call(what=rbind, pbmcapply::pbmclapply(eps_vals, function(x) {
        f(y, x, ...)
      }, mc.cores = n_cores))
    })
  
  stats %>% reshape2::melt(measure.vars=c())
  
}

calc_dic = function(results, eps, use_sl=FALSE, bs=FALSE) {
  
  if (bs) {
    if (use_sl) {
      d = suppressMessages(add_bootstrap_ms_statistics_syn(results, eps=eps, ifelse(isTRUE(bs), 100, bs)))
      d$k = d$k_sl
    } else {
      d = suppressMessages(add_bootstrap_ms_statistics(results, eps=eps, ifelse(isTRUE(bs), 100, bs)))
      d$k = d$k_sl
    }
  } else {
    d = calc_model_stats_syn_lik(results, eps, use_sl=use_sl)
  }
  
  with(d, {
    data.frame(eps=eps, DIC=DIC, pD=pD, pD2=pD2, DIC2=DIC2, k=k)
  })
  
}

calc_aic = function(results, eps, use_sl=FALSE, bs=FALSE) {
  if (bs) {
    if (use_sl) {
      d = suppressMessages(add_bootstrap_ms_statistics_syn(results, eps=eps, ifelse(isTRUE(bs), 100, bs)))
    } else {
      d = suppressMessages(add_bootstrap_ms_statistics(results, eps=eps, ifelse(isTRUE(bs), 100, bs)))
    }
    d$eps = eps 
    d$k = d$k_sl
  } else {
    d = calc_model_stats_syn_lik(results, eps, use_sl=use_sl)
  }
  
  with(d, {
    data.frame(eps=eps, AIC=AIC, AIC2=AIC2, k=k)
  })
}

wrapper_plot_range_of_eps_vs_dic = function(dir, n_cores=1, n_points=50, fname="model_selection_data_v3[.]rds", use_sl=FALSE, ...) {
  data = load_model_selection_data(dir, fname = fname)
  stats = calc_dic_across_eps(data, n_cores, n_points, use_sl, ...)
  plot_dic_values(stats)
}

wrapper_plot_range_of_eps_vs_aic = function(dir, n_cores=1, n_points=50, fname="model_selection_data_v3[.]rds", use_sl=FALSE, ...) {
  data = load_model_selection_data(dir, fname = fname)
  stats = calc_aic_across_eps(data, n_cores, n_points, use_sl, ...)
  plot_stat_values_across_eps(stats, c("AIC","AIC2"))
}



