library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(CHESS)
library(gridExtra)
options(ignore.interactive=TRUE) # always print progress bar.
source("setup_environment.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

stat_used = "AIC2_sl" # statistic to use for model selection
order_of_models = c("Neutral","Selection","Selection + Death", "Selection x 2")
fig_dir = "inference_plots/posteriors_update"
input_trees = readRDS("input_data/final_tree_set.rds")
mutation_rate_cor_file = "input_data/fudge_factors.csv"
excl_cases = c("C543", "C531", " C543", " C531")
dump_data = TRUE

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# overload ggsave function if plot objects should be dumped
if (dump_data) {
  ggsave = function(filename, plot, ...) {
    of = file.path(dirname(filename), paste0(".", gsub("[.][a-zA-Z]+$", ".rds", basename(filename))))
    saveRDS(plot, of)
    ggplot2::ggsave(filename, plot, ...)
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Find model data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

result_sets = 
  file.path(result_dirs, "result_set.rds") %>% # posterior particle files// kind of unused?
  magrittr::set_names(paste0(get_model_set_name(.), " - ", get_model_superset_name(.))) %>% 
  (function(x) x[file.exists(x)])

model_selection_files = # model selection results
  dirname(result_dirs) %>% unique() %>% 
  list.files("model_selection_results_f2[.]rds", recursive=TRUE, full.names=TRUE) %>% 
  magrittr::set_names(paste0(get_model_set_name(.), " - ", get_model_superset_name(.)))

target_trees = # target trees
  dirname(model_selection_files) %>% 
  magrittr::set_names(basename(.)) %>% 
  split(., gsub(" -.*", "", get_model_set_name(.))) %>% 
  lapply(lapply, list.files, pattern="target[.]rds", recursive=TRUE, full.names=TRUE) %>% 
  lapply(lapply, function(x) magrittr::set_names(x, basename(dirname(x)))) %>% 
  lapply(lapply, lapply, readRDS)

n_tips = # number of tips per case
  input_trees %>% 
  lapply(get_n_tip, "cancer") %>% 
  reshape2::melt() %>% 
  magrittr::set_colnames(c("n_tips","case")) %>% 
  dplyr::filter(!case %in% excl_cases)

mutation_rate_fudge_factor = 
  readr::read_csv(mutation_rate_cor_file) %>% 
  as.data.frame() %>% magrittr::set_rownames(.$case)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load datasets ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Loading model selection datasets:\n")

model_selection_data = 
  model_selection_files %>% 
  pbapply::pblapply(load_model_dataset) %>% 
  (function(x) x[sapply(x, NROW) > 0])

wh_stats = c(
  "Delta_AIC2_sl",
  "AIC2_sl",
  "neg_log_lik2_sl",
  "p_model",
  "eps",
  "n_param"
)


cat("Loading model selection posteriors:\n")

.load_dataset = function(f) {
  d = readRDS(f)
  wh = c("particles", "observations", "states")
  d[names(d) %in% wh]
}

result_data = pbapply::pblapply(result_sets, .load_dataset)
names(result_data) = dirname(result_sets)


cat("Calculating statics matrix:\n")

stat_matrix =
  model_selection_data %>% 
  lapply(lapply, "[", c(wh_stats,"dir")) %>%
  lapply(lapply, unlist) %>% 
  lapply(lapply, as.data.frame) %>% 
  lapply(lapply, t) %>% 
  lapply(lapply, data.frame) %>% 
  reshape2::melt(measure.vars = c()) %>% 
  dplyr::mutate(case=sapply(strsplit(L1, " - "), "[[", 2)) %>% 
  dplyr::mutate(set=gsub(" - .*", "", L1)) %>% 
  dplyr::mutate(superset=get_model_superset_name(dir)) %>% 
  dplyr::filter(!case %in% excl_cases) %>% 
  dplyr::mutate(model = 
    factor(L2, unique(L2[order(n_param)]), gsub(" [+] Push", "", format_model_string(unique(L2[order(n_param)]), FALSE)
  ))) %>%
  merge(n_tips, all.x=TRUE)
  
for (stat in wh_stats) {
  stat_matrix[,stat] = as.numeric(stat_matrix[,stat])
}

stat_matrix = # add delta AIC data
  stat_matrix %>% 
  split(., paste0(.$superset, " - ", .$set, " - ", .$case)) %>%
  lapply(function(x) {
    x$AIC2_sl = as.numeric(x$AIC2_sl)
    x$Delta_AIC2_sl = x$AIC2_sl - min(x$AIC2_sl)
    x$Delta_AIC2_sl_undecided = sum(x$Delta_AIC2_sl < 4) > 1
    return(x)
  }) %>% do.call(what=rbind)


get_best_model_line = function(x, stat) {
  x %>% 
    dplyr::mutate(model=as.character(model)) %>%
    split(., paste0(.$case, "-", .$set)) %>%
    lapply(arrange, as.numeric(get(stat))) %>%
    lapply(function(x) {
      x = head(x, 1)
      ud_col = paste0(stat, "_undecided")
      if (ud_col %in% colnames(x)) {
        if (x[,ud_col]) x$model = "Undecided"
      }
      return(x)
    }) %>% 
    do.call(what=rbind)
}

saveRDS(model_selection_data, file.path(fig_dir, "ms_data.rds"))
saveRDS(stat_matrix, file.path(fig_dir, "stat_matrix.rds"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Export summary tables
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


begin{table}[ht]
\centering
\small
\caption{Model selection Variable Sampling - WGS + LP.} 
\begin{tabular}{rlllrlrl}
\hline
& \textbf{Case} & \textbf{Best Model} & \textbf{$\Delta$AIC-N} & \textbf{$\Delta$AIC-S} & \textbf{$\Delta$AIC-Sx2}& \textbf{p} & \textbf{Comment}\\ 
\hline
1 & C516 & \textbf{Selection} & 2.41 & 0.00 & 7.57 & 0.25 &  Selection in B (Unknown driver?) \\ 
2 & C518 & \textbf{Selection x 2} & 3.49 & 1.31 & 0 & 0.03 & Selection in A \& B (PTEN p.C136R) \\ 
5 & C524 & \textbf{Selection} & 2.1 & 0.00 & 4.79 & 0.05 & Selection in B (PIK3CA p.C378R) \\ 
6 & C525 & \textbf{Selection} & \textbf{9.01} & 0.00 & 5.91 & 0.03 & Selection in C (PIK3CA p.Q546P) \\ 
8 & C528 & Neutral & 0 & 2.60 & - & 0.00 &  \\ 
9 & C530 & Neutral & 0 & 5.05 & - & 0.01 &  \\ 
10 & C531 & \textbf{Selection} & 3.74 & 0.00 & 6.51 & 0.03 & Selection in D (PIK3CA p.Q546K) \\ 
11 & C532 & Neutral & 0 & 10.30 & 20.6 & 0.13 &  \\ 
12 & C536 & Neutral & 0 & 15.90 & - & 0.59 &  \\ 
13 & C537 & Neutral & 0 & 6.19 & - & 0.01 &  \\ 
14 & C538 & \textbf{Selection} & \textbf{14.6} & 0.00 & 6.29 & 0.07 & Selection in D (RNF43 p.Q153*) \\ 
15 & C539 & \textbf{Selection} & \textbf{28.6} & 0.00 & 3.05 & 0.03 & Selection in A (KRAS p.G12C) \\ 
16 & C542 & \textbf{Selection} & \textbf{18.4} & 0.00 & 1.56 & 0.00 & Selection in A,B \& D? (chr1p loss?) \\ 
17 & C543 & Neutral & 0 & 0.88 & - & 0.00 &  \\ 
18 & C544 & Neutral & 0 & 9.77 & - & 0.22 &  \\ 
20 & C548 & Neutral  & 0 & 7.28 & - & 0.19 &  \\ 
21 & C549 & \textbf{Selection} & 3.16 & 0.00 & 7.06 & 0.03 & Selection in A(chr1p loss?)  \\ 
22 & C550 & Neutral & 0 & 8.58 & - & 0.03 &  \\ 
23 & C551 & \textbf{Selection} & 1.59 & 0.00 & 7.52 & 0.09 & Selection in A \& B (Unknown driver?) \\ 
24 & C552 & Neutral & 0 & 11.30 & - & 0.50 &  \\ 
25 & C554 & Neutral & 0 & 4.08 & - & 0.12 &  \\ 
26 & C555 & Neutral & 0 & 8.16 & 16.3 & 0.55 &  \\ 
27 & C559 & \textbf{Selection} & \textbf{5.05} & 0.00 & - & 0.17 & Selection  in A \& B (Unknown driver?) \\ 
28 & C560 & Neutral & 0 & 10.50 & - & 0.07 &  \\ 
29 & C561 & Neutral & 0 & 5.09 & - & 0.19 &  \\ 
30 & C562 & Selection & 0.321 & 0.00 & - & 0.12 & \\ 
31 & C531 (excl. D1-G7) & \textbf{Selection} & \textbf{8.9} & 0.00 & 5.02 & 0.00 & Selection in B (SMAD4 p.A118V) \\ 
32 & C543 (excl. A1-G9) & Neutral & 0 & 4.17 & 11.6 & 0.01 &  \\ 
\hline
\end{tabular}
\end{table}


get_table = function(x) {
  # get delta aic table
  delta_aic_tab = 
    with(matrix_view, tapply(Delta_AIC2_sl, list(case, model), c)) %>% 
    as.data.frame() %>% 
    magrittr::set_colnames(paste0("Delta_", colnames(.)))
  
  # get best model
  best_model = split(matrix_view, matrix_view$case) %>% 
    lapply(get_best_model_line, "AIC2_sl") %>% 
    do.call(what=rbind) %>% 
    dplyr::mutate(i=seq_along(case)) %>% 
    dplyr::mutate(comment="") %>% 
    dplyr::mutate(p_model = signif(as.numeric(p_model, 3)))
  
  # create table
  cbind(best_model, delta_aic_tab[as.character(best_model$case),]) %>% 
    dplyr::select(case, model, neg_log_lik2_sl, AIC2_sl,  Delta_Neutral, Delta_Selection, `Delta_Selection x 2`, p_model, comment) %>% 
    magrittr::set_colnames(c("Case","Best Model", "NLL", "AIC", "$\\Delta$AIC-N", "$\\Delta$AIC-S", "$\\Delta$AIC-Sx2", "P", "Comment"))
}


ss =  paste0(stat_matrix$superset, " - ", stat_matrix$set)
for (ss_ in unique(ss)) {
  # tabulate data
  tab = get_table(stat_matrix[ss == ss_,])
  tab$Case = gsub("[.].*", "", tab$Case)
  rownames(tab) = NULL
      
  # highlight non-neutral delta aic >4
  wh = which(as.numeric(tab[,"$\\Delta$AIC-N"]) > 4)
  tab[,"$\\Delta$AIC-N"] = round(tab[,"$\\Delta$AIC-N"], 3)
  tab[wh,"$\\Delta$AIC-N"] = paste0("BOLD(", tab[wh,"$\\Delta$AIC-N"], ")")

  # write table
  bold = function(x) gsub('BOLD[(](.*)[)]',paste0('\\\\textbf{\\1','}'),x)
  tab_tex = xtable::xtable(tab, paste0("Model selection ", gsub("&", "+", ss_), "."))
  ofile = file.path(fig_dir, "tex_tables", paste0(gsub("[ &]", "_", ss_), ".tex"))
  dir.create(dirname(ofile), FALSE, TRUE)
  print(tab_tex, file = ofile, compress = FALSE, sanitize.text.function = bold)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Summary of number of tips per model ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

for (stat in c("AIC2_sl", "Delta_AIC2_sl")) {
  
  for (ss in unique(stat_matrix$superset)) {
    
    best_models =
      stat_matrix %>% 
      dplyr::filter(superset == ss) %>%
      get_best_model_line(stat) %>%
      dplyr::select(model, n_tips, set)
    
    
    bootstrap_p = function(x, y) {
      set.seed(123)
      gr = c(rep(1, length(x)), rep(2, length(y)))
      v = c(x, y)
      delta = (diff(tapply(v, gr, mean)[c("1","2")]))
      delta_bs = replicate(10000, (diff(tapply(v, sample(gr), mean)[c("1","2")])))
      list(p.value=mean((delta) > (delta_bs)))
    }
    
    plot =
      best_models %>%
      dplyr::mutate(model=gsub(" .*", "", model)) %>%
      dplyr::mutate(set = paste0(gsub(" only", "", set), " tips*")) %>%
      dplyr::mutate(set = factor(set, c("WGS tips*","WGS & LP tips*"))) %>%
      ggplot(aes(x=model, y=n_tips)) +
        geom_boxplot() +
        facet_wrap(~set, scales = "free_x") +
        geom_jitter(width=0.1, height=0) +
        ggsignif::stat_signif(comparisons = list(c("Selection","Neutral")), test="bootstrap_p") +
        xlab("Classification") +
        ylab("Number of tips") +
        scale_y_continuous(expand = expand_scale(mult=c(0, 0.15)), n.breaks = 4, limits=c(0, NA)) +
        labs(caption = "*Cancer samples only") +
        theme(plot.caption = element_text(size=9))

    outfile = file.path(fig_dir, ss, paste0("classification_vs_tips_", ss, "_", stat, ".pdf"))
    dir.create(dirname(outfile), FALSE, TRUE)
    ggsave(outfile, plot, width=ifelse(grepl("Delta", stat), 3, 2.5), height=2.7)

  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Summary of number of selected models ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

for (ss in unique(stat_matrix$superset)) {
  
  # get best models for various stats
  stats = c("AIC"="AIC2_sl","Delta~AIC"="Delta_AIC2_sl")
  best_models = NULL
  for (i in seq_along(stats)) {
    
    best_models_ =
      stat_matrix %>% 
      dplyr::filter(superset == ss) %>%
      get_best_model_line(stats[i]) %>%
      dplyr::select(model, set, case) %>% 
      dplyr::mutate(stat=names(stats)[i]) %>% 
      dplyr::mutate(model=as.character(model))
    
    best_models = 
      rbind(best_models, best_models_)
    
  }

  # plot model numbers
  plot_model_selection_number = 
    best_models %>% 
    dplyr::filter(!case %in% excl_cases) %>% 
    dplyr::mutate(model = factor(model, c("Undecided", rev(order_of_models)))) %>% 
    dplyr::mutate(set=factor(set, c("WGS only","WGS & LP"))) %>% 
    dplyr::mutate(stat=factor(stat, names(stats))) %>% 
    ggplot(aes(x=stat, fill=model)) + 
    geom_bar() + 
    facet_wrap(~set, scales = "free") + 
    xlab("Statistic") + 
    ylab("Cases") + 
    labs(fill="Model") + 
    scale_x_discrete(
      labels=parse(text=c("AIC","Delta~AIC")),
      breaks=c("AIC","Delta~AIC")
    ) + 
    scale_fill_manual(
      values = rev(c("gray60", scales::brewer_pal("seq", "Set2", 1)(3))), 
      breaks = rev(c("Undecided", order_of_models[-3]))
    ) +
    ggtitle(ss)
  
  ofile = paste0(file.path(fig_dir, ss, "model_selection_total_numbers_aic_delta_aic.pdf"))
  dir.create(dirname(ofile), FALSE, TRUE)
  ggsave(ofile, plot_model_selection_number, width=4, height=2.8)
  
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plots of posterior using best model ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

get_tree_height = function(x) {
  ids = x$tip.label[x$tip.label != "GL"]
  a = THmisc::annotation_from_barcode(ids, TRUE)
  x = ape::keep.tip(x, c("GL", ids[a$tissue_type=="cancer" & a$analyte_name == "WGS"]))
  d = cophenetic(x)
  n_clonal = min(d["GL", d["GL", ] > 0])
  mean(d["GL", d["GL", ] > 0] - n_clonal)
}

scale_m_rates = function(d, ff, max_mut=250) {
  
  if ("map" %in% colnames(d)) {
    d %<>%
      dplyr::mutate(map=if_else("mutation_rate" == parameter, map/ff[gsub("([.].*| )", "", case)], map)) %>% 
      dplyr::mutate(map=if_else("mutation_rate" == parameter, scales::squish(map, c(0, max_mut)), map)) %>% 
      dplyr::mutate(shape=if_else("mutation_rate" ==  parameter & map >= max_mut, 17, 19))
  }
  
  if ("lower" %in% colnames(d)) {
    d %<>%
      dplyr::mutate(lower=if_else("mutation_rate" ==  parameter, lower/ff[gsub("([.].*| )", "", case)], lower)) %>% 
      dplyr::mutate(lower=if_else("mutation_rate" ==  parameter,  scales::squish(lower, c(0, max_mut)), lower))
  }
  
  if ("upper" %in% colnames(d)) {
    d %<>%
      dplyr::mutate(upper=if_else("mutation_rate" == parameter, upper/ff[gsub("([.].*| )", "", case)], upper)) %>% 
      dplyr::mutate(upper=if_else("mutation_rate" == parameter, scales::squish(upper, c(0, max_mut)), upper))
  }
  
  
  if ("value" %in% colnames(d)) {
    d %<>%
      dplyr::mutate(value=if_else("mutation_rate" == parameter, value/ff[gsub("([.].*| )", "", case)], value)) %>%
      dplyr::mutate(value=if_else("mutation_rate" == parameter, scales::squish(value, c(0, max_mut)), value))
  }
  
  return(d) 
}

get_case_order = function(map) {
  
  st2 = tryCatch(
    map["clone_start_times.2", ],
    error = function(e)
      rep(1, NCOL(map))
  )
  
  st3 = tryCatch(
    map["clone_start_times.3", ],
    error = function(e)
      rep(1, NCOL(map))
  )
  
  pp = map["push_power.1", ]
  
  case_order =  colnames(map)[order(st3, st2, pp > 5, pp, na.last = FALSE)] #unique(pp$case[order(pp$model, pp$map>)])
  case_order = case_order[grepl("[ ]", case_order)]
  
  c(gsub(" ", "", case_order), case_order)
}

add_pp_insert_data_fit = function(d, pp_low_freq) {
  
  new_level = c(levels(d$parameter_n), "d[push] ")

  fits_of_best_model_pp = d %>% 
    dplyr::filter(parameter == "push_power.1") %>% 
    dplyr::mutate(parameter = "push_power.1_low_freq") %>% 
    dplyr::mutate(parameter_n = factor("d[push] ", new_level)) %>%
    dplyr::mutate(shape = if_else(map >= pp_low_freq, 17, 19)) %>% 
    dplyr::mutate(map = scales::squish(map, c(0, pp_low_freq))) %>% 
    dplyr::mutate(upper = scales::squish(upper, c(0, pp_low_freq))) %>% 
    dplyr::mutate(lower = scales::squish(lower, c(0, pp_low_freq)))
  
  d %>% 
    dplyr::mutate(parameter_n = factor(parameter_n, new_level)) %>% 
    rbind(fits_of_best_model_pp)
  
}

add_pp_insert_data_post = function(d, pp_low_freq) {
  
  new_level = c(levels(d$parameter_n), "d[push] ")

  post_dist_data_low_pp = d %>% 
    dplyr::filter(parameter == "push_power.1") %>% 
    dplyr::filter(value <= pp_low_freq) %>% 
    dplyr::mutate(parameter = "push_power.1_low_freq") %>% 
    dplyr::mutate(parameter_n = factor("d[push] ", new_level))
    
  d %>% 
    dplyr::mutate(parameter_n = factor(parameter_n, new_level)) %>% 
    rbind(post_dist_data_low_pp)

}

plot_posterior = function(fit, dist_data=NULL) {
  
  plot_model_fits = 
    fit %>% 
    dplyr::filter(!is.na(model)) %>% 
    dplyr::filter(parameter != "region_scale_factor.1") %>%
    ggplot(aes(x=case))
  
  if (!is.null(dist_data)) {
    plot_model_fits = plot_model_fits + 
      geom_violin(
        data = dist_data,
        aes(weight = weight, y = value, fill = model),
        alpha = 1,
        scale = "width",
        width = 0.8
      )
  }
  
  plot_model_fits = plot_model_fits + 
    geom_point(aes(y=map,  color=model, shape=shape), color="gray10") + 
    geom_linerange(aes(y=map, ymin=lower, ymax=upper, color=model), color="gray10", alpha=0.8) + 
    facet_grid(parameter_n~set, scales = "free", space = "free_x", labeller = labeller(parameter_n=label_parsed)) + 
    xlab("") + 
    ylab("Posterior") + 
    labs(alpha="Good fit", shape="B2F estimates", color="Selected model", fill="Selected model") + 
    scale_color_grey() + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
    theme(strip.text.y = element_text(angle=0, hjust=0.5, vjust=0.5)) + 
    theme(legend.position="bottom") + 
    theme(legend.box = "vertical") + 
    background_grid() + 
    scale_fill_brewer(palette = "Set2") + 
    scale_shape_identity() + 
    theme(strip.background = element_rect(fill="gray10"), strip.text = element_text(color="gray90"))
  
  return(plot_model_fits)
}

# variables
max_mut = 500
pp_low_freq = 50
tree_heights = sapply(input_trees, get_tree_height)
ff_lookup_table = with(mutation_rate_fudge_factor, magrittr::set_names(fudge_factor, case))


for (stat in c("AIC2_sl", "Delta_AIC2_sl")) {
  for (ss in unique(stat_matrix$superset)) {
  
    # get summary stats of fits
    fits_of_best_model =
      lapply(model_selection_data, get_best_model_fit, statistic = stat) %>%
      reshape::melt(measure.vars = c()) %>%
      dplyr::mutate(case = sapply(strsplit(L1, " - "), "[[", 2)) %>%
      dplyr::mutate(set = gsub(" - .*", "", L1)) %>%
      scale_m_rates(ff_lookup_table, max_mut) %>%
      dplyr::mutate(case = ifelse(set == "WGS & LP", paste0(" ", case), case)) %>% # prefix
      dplyr::mutate(model = factor(format_model_string(model))) %>%
      dplyr::mutate(parameter_n = factor(params[parameter], params)) %>%
      dplyr::mutate(superset = get_model_superset_name(dir)) %>%
      dplyr::mutate(set = factor(set, c("WGS only", "WGS & LP"))) %>%
      dplyr::filter(!is.na(parameter_n)) %>%
      dplyr::filter(superset == ss) %>%
      dplyr::filter(parameter != "region_scale_factor.1") %>% 
      add_pp_insert_data_fit(pp_low_freq)
  
    # order models
    model_order = unique(order_of_models, as.character(fits_of_best_model$model))
    fits_of_best_model$model = factor(fits_of_best_model$model, model_order, ordered = TRUE)
  
    # order cases
    map_vals = tapply(fits_of_best_model$map, list(fits_of_best_model$parameter, fits_of_best_model$case), mean)
    case_order = get_case_order(map_vals)
    fits_of_best_model$case = factor(fits_of_best_model$case, case_order, ordered = TRUE)
  
    # get full posterior distribution
    post_dist_data = get_post_dist_data(fits_of_best_model)
    if (!is.null(post_dist_data)) {
    post_dist_data = post_dist_data %>% 
      scale_m_rates(ff_lookup_table, max_mut) %>%
      dplyr::filter(parameter != "region_scale_factor.1") %>% 
      dplyr::mutate(model = factor(model, model_order, ordered = TRUE)) %>% 
      dplyr::mutate(case = factor(case, case_order, ordered = TRUE)) %>%
      add_pp_insert_data_post(pp_low_freq)
    }
  
    # plot the results
    plot_model_fits = 
      plot_posterior(
        fit = fits_of_best_model, 
        dist_data = post_dist_data
      )
    
    ofile = paste0(file.path(fig_dir, ss, "posterior_model_fits", paste0("posterior_model_fits", stat, ".pdf")))
    dir.create(dirname(ofile), FALSE, TRUE)
    ggsave(ofile, plot_model_fits, width=9, height=8)
    
  }
}


# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# # Plots of selection projected on trees ####
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# # 
# # 
# for (cset in unique(plot_model_fits$data$set)) {
#   try({
#     
#     # create panel of trees in order of model and bd growths:
#     pp_order = 
#       plot_model_fits$data %>% 
#       dplyr::filter(parameter == "push_power.1") %>% 
#       dplyr::arrange(-as.numeric(model), -map) %>% 
#       dplyr::filter(!gsub(" ", "", case) %in% excl_trees)
#     
#     case_ord = 
#       pp_order %>% 
#       dplyr::select(case) %>% 
#       unlist() %>% 
#       as.character()
#     
#     post_particle_files = 
#       file.path(pp_order$dir, paste0("posterior_particles_v", .posterior_particle_version, ".rds")) %>% 
#       magrittr::set_names(pp_order$case)
#     
#     
#     # create plots of trees:
#     trees = input_trees[gsub(" ", "", case_ord)] %>% magrittr::set_names(case_ord)
#     if (cset == "WGS only")  trees = lapply(trees, function(x) ape::drop.tip(x, x$tip.label[grepl("Added", x$tip.label)]))
#     #trees = trees %>% lapply(remove_root_tip, "GL")
#     
#     cat("Creating tree plots...\n")
#     
#     tree_plots = 
#       pbmcapply::pbmclapply(case_ord, function(case) {
#         
#         plt = 
#           plot_tree_with_edge_freq(
#             dfile = post_particle_files[[case]], 
#             tree_use = trees[[case]], 
#             size_lty=1
#           )
#         
#         # add driver labels?
#         if (gsub(" ", "", case) %in% names(sc_driver_annotation)) {
#           plt = 
#             add_driver_labels_to_tree(
#               plt, 
#               sc_driver_annotation[[gsub(" ", "", case)]]
#             )
#         }
#         
#         return(plt)
#       }, mc.cores = 1)
#     
#     names(tree_plots) = case_ord
#     legend = cowplot::get_legend(tree_plots[[1]]) 
#     
#     # modify labels of trees
#     for (i in seq_along(tree_plots)) {
#       
#       params_i = 
#         plot_model_fits$data %>% 
#         dplyr::filter(case == names(tree_plots)[i])
#       
#       pp = signif(dplyr::filter(params_i, parameter == "push_power.1")$map, 2)
#       model = as.character(dplyr::filter(params_i, parameter == "push_power.1")$model)
#       max_x = max(tree_plots[[i]]$data$x)
#       pp_color = c("1"="#339900","2"="#99cc33","3"="#f66630","4"="#cc3300")[as.character(as.numeric(cut(pp, c(0, 5, 10, 50, Inf))))]
#       model_color = c("Selection x 2" ="#cc3300", "Selection"="#f04c35", "Neutral"="#339900")[model]
#       
#       data_pp = data.frame(x=max_x, y=-0.5, pp=paste0("d_push = ", pp))
#       
#       tree_plots[[i]] = 
#         tree_plots[[i]] + 
#         geom_text(data=data_pp, aes(x=x*0.8, y=y, label=pp), size=3, inherit.aes=FALSE, hjust=0, color=pp_color) + 
#         ggtitle(paste0(names(tree_plots)[i], " (", model, ")")) + 
#         theme(plot.title = element_text(color = model_color)) + 
#         guides(colour_new=FALSE)
#     } 
#     
#     plots = arrangeGrob(grobs=do.call(grid::gList, lapply(c(tree_plots, list(ggpubr::as_ggplot(legend))), ggplotGrob)), nrow=6, ncol=5)
#     ofile = file.path(fig_dir, ss, "sorted_trees", paste0(gsub("[ ]", "", gsub("&", "", cset)), "_", stat, ".pdf"))
#     dir.create(dirname(ofile), FALSE, TRUE)
#     ggsave(ofile, plots, width=10, height=12)
#     
#     
#   })
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plots of posterior using best model ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 

for (stat in c("AIC2_sl", "Delta_AIC2_sl")) {
  for (ss in unique(stat_matrix$superset)) {
    
    msi_positiv = c("C536", "C548", "C516", "C518", "C552", "C562")
    
    # get summary stats of fits
    mrate =
      lapply(model_selection_data, get_best_model_fit, statistic = stat) %>%
      reshape::melt(measure.vars = c()) %>%
      dplyr::mutate(case = sapply(strsplit(L1, " - "), "[[", 2)) %>%
      dplyr::mutate(set = gsub(" - .*", "", L1)) %>%
      scale_m_rates(ff_lookup_table, Inf) %>%
      dplyr::mutate(case = ifelse(set == "WGS & LP", paste0(" ", case), case)) %>% # prefix
      dplyr::mutate(model = factor(format_model_string(model))) %>% 
      dplyr::filter(parameter == "mutation_rate") %>% 
      dplyr::filter(get_model_superset_name(dir) == ss) %>%
      split(paste0(.$case, .$set)) %>% 
      lapply(function(x){stopifnot(NROW(x)==1); x}) %>% 
      do.call(what=rbind) %>% 
      group_by(gsub("( |[.].*)", "", case) %in% msi_positiv, set) %>% 
      dplyr::summarise(median_m = median(map, na.rm=TRUE), mean_m = mean(map, na.rm=TRUE))
    
    print(paste0(stat, " - ", ss))
    print(mrate)
    
  }
}
