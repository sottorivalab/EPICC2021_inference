
get_model_set_name = function(x) {
  
  sn = case_when(
    grepl("excluding_lp_random_pushing", x) ~ "WGS only (old)",
    grepl("wgs_only_pushing_to_edge", x) ~ "WGS only",
    grepl("including_lp_pushing_to_edge", x) ~ "WGS & LP"
  )
  
  paste0(sn, " - ", basename(dirname(x)))
}

get_model_superset_name = function(x) {
  
  sn = case_when(
    grepl("alma_5000x200x1_alpha_0.95_overdispersed_v3_popsize_100000", x) ~ "Variable Sampling",
    grepl("alma_500x100x25_alpha_0.95/", x) ~ "Fixed Sampling",
    grepl("alma_500x100x25_alpha_0.95_overdispersed/", x) ~ "Fixed Sampling - Overdispersed",
    grepl("alma_500x100x25_alpha_0.95_clip_tips/", x) ~ "Fixed Sampling - Cliped tips",
    grepl("including_lp_pushing_to_edge", x) ~ "Other"
  )
  
  paste0(sn)
}

get_model_fit = 
  function(x, model) {
    
    if (!model %in% names(x)) return(NULL)
    
    # get summary stats:
    estimated_summaries = c("map","lower","upper")
    est = x[[model]][estimated_summaries]
    params = names(est[[1]])
    stopifnot(all(sapply(est, function(x) all(names(x) == params))))
    est_merged = cbind(data.frame(parameter=params, do.call(cbind, est)))
    
    return(est_merged)
  }

get_best_model_fit = 
  function(x, statistic="DIC") {
    
    # find best model:
    stat = sapply(x, "[[", statistic)
    best = which.min(stat)
    
    # add model label to data
    est_merged = get_model_fit(x, names(stat)[best])
    est_merged$model = names(stat)[best]
    est_merged$dir = x[[best]]$dir
    
    return(est_merged)
  }

load_model_dataset = function(x) {
  
  els = c("d_dhat", "map_ndim", "map", "mean", "median", "lower", "upper", 
          "est_dhat", "particle_dhat", "d_dhat_min_sc_frac", "d_dhat_permutated", 
          "particles_dbar", "d_dbar", "d_dbar_min_sc_frac", "particles_dbar2", 
          "d_dbar2", "d_dbar2_min_sc_frac", "n_param", "neg_log_lik", "AIC", 
          "neg_log_lik2", "AIC2", "Dbar", "Dhat", "pD", "pV", "DIC", "Dbar2", 
          "pD2", "pV2", "DIC2", "p_model", "M", "N_p", "N", "N_sim", "eps", 
          "eps_min_smc", "eps_noise", "d_dhat_syn", "d_dbar_syn", "d_dbar2_syn", 
          "neg_log_lik_sl", "neg_log_lik2_sl", "neg_log_lik_map_sl", "k_sl", 
          "Dbar_sl", "Dbar2_sl", "Dhat_sl", "pV_sl", "pV2_sl", "AIC_sl", 
          "AIC2_sl", "pD_sl", "pD2_sl", "DIC_sl", "DIC2_sl",
          "Delta_AIC_sl","Delta_AIC2_sl","dir")
  
  d = readRDS(x)
  if (length(d) == 0) return(NULL) 
  
  for (i in seq_along(d)) {
    d[[i]][duplicated(names(d[[i]]))]=NULL
    d$kernel = NULL
    d[[i]]$dir = file.path(dirname(x), names(d)[i])
  }
  
  for (stat in c("AIC_sl", "AIC2_sl")) {
    min_stat = min(sapply(d, "[[", stat), na.rm=TRUE)
    k = which.min(sapply(d, "[[", "n_param"))
    
    for (i in seq_along(d)) {
      if (i %in% k) {
        d[[i]][[paste0("Delta_", stat)]] = d[[i]][[stat]]
      } else {
        d[[i]][[paste0("Delta_", stat)]] = d[[i]][[stat]] + 4
      }
    }
  }
  
  
  for (i in seq_along(d)) {
    d[[i]] = d[[i]][els]
  }
  
  d
}

get_post_dist_data = function(x) {
  
  if (!exists("result_data")) return(NULL)
  res_data  = NULL
  sets = unique(x$dir)
  par = unique(x$parameter)
  
  for (set in sets){
    
    if (!set %in% names(result_data)) next()
    eps = min(result_data[[set]]$states$epsilon)
    d_cur = recalculate_fracs_below_eps(result_data[[set]], eps)$particles
    
    par_cur = c(par[par %in% colnames(d_cur)], "frac_below_eps")
    d_cur = d_cur[,par_cur]
    d_cur$frac_below_eps = d_cur$frac_below_eps / sum(d_cur$frac_below_eps, na.rm=TRUE)
    
    wh = !is.na(d_cur$frac_below_eps)
    d_cur_l = reshape2::melt(d_cur[wh,], id.vars=c("frac_below_eps"))
    colnames(d_cur_l) = c("weight","parameter","value")
    
    other_vals = as.list(head(x[x$dir == set,], n=1))
    d_cur_l$parameter_n = factor(params[as.character(d_cur_l$parameter)], params)
    d_cur_l$case = other_vals[["case"]]
    d_cur_l$set = other_vals[["set"]]
    d_cur_l$model = other_vals[["model"]]
    
    res_data = rbind(res_data, d_cur_l)
  }
  
  res_data
}

get_n_tip = function(x, tissue=NULL) {
  lab = x$tip.label[x$tip.label!="GL"]
  if (!is.null(tissue)) {
    annot = suppressWarnings(THmisc::annotation_from_barcode(lab, TRUE))
    lab = lab[as.character(annot$tissue_type) %in% tissue]
  }
  length(lab)
}
