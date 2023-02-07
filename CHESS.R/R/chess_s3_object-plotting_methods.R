#' Plot population sizes over time for a CHESS_S3 simulation.
#'
#' @param x A CHESS_S3 object.
#'
#' @return A ggplot object.
#' @export
#' 
#' @examples plot_history_CHESS_S3(new_simulation(20)) # a small 20x20x1 simulation
plot_history_CHESS_S3 = function(x) {
  
  clone_colors = c(
    "1" = "#377EB7",
    "2" = "#E3211C",
    "3" = "#4DAE4A",
    "4" = "#974EA2",
    "5" = "#F680BE",
    "6" = "#FE7F00"
  )
  
  checkmate::assertClass(x, "CHESS_S3")
  
  x = add_history_CHESS_S3(x)
  
  history = x$history
  colnames(history) = gsub("type", "", colnames(history))
  
  # point of introduction
  idx_start = apply(history[,-1,drop=FALSE] > 0, 2, function(x) {wh=which(x); if (length(wh)) return(wh[1]) else return(NROW(history)) })
  start_time = history$time[idx_start]
  start_data = data.frame(time=start_time, variable=colnames(history)[-1])
  
  
  # subset to fewer values:
  count_data = history[,-1,drop=FALSE]
  n_points = NROW(history)
  keep_low_count = apply(count_data > 0 & count_data < 100, 1, any)     # keep points with low counts (>0 & <100)
  keep_n_element = seq_along(keep_low_count) %% (n_points %/% 200) == 0 # keep 200 points in general
  history = history[keep_low_count | keep_n_element, ]
  
  
  # time line
  total_data = data.frame(time=history[,1], value=apply(history[,-1, drop=FALSE], 1, sum))
  time_line_data = history %>% reshape::melt(id.vars="time")
  
  relative_data = cbind(time=history[,1], history[,-1,drop=FALSE] / total_data$value)
  colnames(relative_data) = colnames(history)
  relative_data = relative_data %>% reshape::melt(id.vars="time")
  
  
  # create the plot
  if (NROW(time_line_data) == 1) {
    
    # return empty plot if data are missing
    time_line_plot = 
      ggplot(data=data.frame()) + 
      geom_text(aes(label="Missing history data", x=0.5, y=0.5)) + 
      xlab("Time") + 
      ylab("Population size") + 
      theme(panel.background = element_rect(fill="transparent",colour = NA)) + 
      theme(plot.background = element_rect(fill="transparent",colour = NA)) + 
      theme(axis.text = element_blank(), axis.ticks = element_blank())
    
  } else {
    
    time_line_plot =
      time_line_data %>%
      ggplot(aes(x = time, y = value, color = variable)) +
      cowplot::theme_cowplot() +
      geom_line(
        data = total_data,
        aes(group = 1),
        color = 'black',
        linetype = 3,
        alpha = 0.9
      ) +
      geom_point(
        data = start_data,
        y = max(time_line_data$value) * -0.03,
        shape = 17
      ) +
      geom_line() +
      xlab("Gillespie time") +
      ylab("Population size") +
      labs(color = "Clone") +
      scale_y_log10() +
      theme(panel.background = element_rect(fill = "transparent", colour = NA)) +
      theme(plot.background = element_rect(fill = "transparent", colour = NA)) +
      theme(legend.position = "bottom")
   
    
    if (!all(na.omit(time_line_data$variable) %in% names(clone_colors))) {
        
        # fall-back option
        time_line_plot = 
          time_line_plot + 
          scale_color_brewer(
            palette = "Set1",
            direction = -1,
            na.value = "white",
            na.translate = FALSE
          )
        
      } else {
        
        time_line_plot = 
          time_line_plot + 
          scale_color_manual(
            breaks = names(clone_colors),
            values = clone_colors,
            na.value = "white",
            na.translate = FALSE
          )
        
      }
    
  }
  
  return(time_line_plot)
  
}

#' Plot space of a CHESS_S3 simulation.
#'
#' @param x  A CHESS_S3 object.
#' @param alpha_burden (optional) flag indicating if the alpha value should indicate the mutation burden of cells (default: true).
#' @param mm_per_dot (optional) value that scales the grid points to a pseudo mm value (default: 0.1).
#' @param labeller_function (optional) function that modifies the tip.labels (default: NULL).
#' @param rotate_by (optional) numeric value indicating by how many radians the space should be rotated (default: 0).
#' @param mirror (optional) flag indicating if the space should be mirrored along the x axis (default: false).
#' @param add_transformation_locations 
#' @param ... 
#'
#' @return A ggplot object.
#' @export
#'
#' @examples sim = new_simulation(20, seed=1) # a small 20x20x1 simulation.
#' @examples sim_plus_sample = add_all_samples(sim, sample_list="random_cells")
#' @examples plot_space_CHESS_S3(sim) 
#' @examples plot_space_CHESS_S3(sim, rotate_by=0.5*pi) # as above, but rotated 90° to the left.
#' @examples plot_space_CHESS_S3(sim, mirror=TRUE) # as above, but mirrored along the x-axis.
#' @examples plot_space_CHESS_S3(sim_plus_sample) # plot with samples
#' @examples labeller_function = function(x) LETTERS[seq_along(x)]
#' @examples plot_space_CHESS_S3(sim_plus_sample, labeller_function=labeller_function) # plot with modified labels
plot_space_CHESS_S3 = function(x, alpha_burden=TRUE, mm_per_dot=0.1, labeller_function=NULL, rotate_by=0, mirror=FALSE, add_transformation_locations=TRUE,  recenter=FALSE, add_polygons=FALSE, ...) {
  
  clone_colors = c(
    "1" = "#377EB7",
    "2" = "#E3211C",
    "3" = "#4DAE4A",
    "4" = "#974EA2",
    "5" = "#F680BE",
    "6" = "#FE7F00"
  )
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Internal helper function ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  sample_coords_to_polygon = function(d) {
    
    if ("sample_coords" %in% names(d)) coords = d[["sample_coords"]] else coords = d
    coords[4:6] = coords[4:6] + 1
    
    xv = coords[c(1,1,4,4,1)] 
    yv = coords[c(2,5,5,2,2)]
    xlab = max(coords[c(1,4)]) #mean(coords[c(1,4)]) # annotation position
    ylab = mean(coords[c(2,5)])
    data.frame(x=xv, y=yv, xlab=xlab, ylab=ylab)
  }
  
  add_location_polygon = function(plot, center, diameter, label) {
    
    if (length(center) == 2) center = c(center, 1) # add missing z coordinate
    if (length(diameter) == 2) diameter = c(diameter, 1) # add missing z coordinate
    
    coords = ceiling(c(center - diameter/2, center + diameter/2 -1))
    polygon_paths = sample_coords_to_polygon(coords)
    polygon_paths$L1 = label
    
    sample_labels = 
      polygon_paths %>% 
      dplyr::select(L1, xlab, ylab) %>% 
      unique()
    
    space_plot = 
      plot + 
      geom_path(aes(x=x,y=y, group=L1), data=polygon_paths, inherit.aes=0) +
      geom_text_repel(aes(label=L1, x=xlab, y=ylab), 
                      data=sample_labels, inherit.aes=0, 
                      min.segment.length=0, force = 1)
    
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assertFlag(alpha_burden, add=al)
  checkmate::assertNumeric(mm_per_dot, any.missing = FALSE, null.ok = FALSE, finite = TRUE, len=1, lower=0, add=al)
  checkmate::assert_true(mm_per_dot > 0.0, add=al)
  checkmate::assert_function(labeller_function, null.ok = TRUE, add=al)
  checkmate::assertNumeric(rotate_by, any.missing = FALSE, null.ok = FALSE, finite = TRUE, len=1, add=al)
  checkmate::assertFlag(mirror, add=al)
  checkmate::assertFlag(add_transformation_locations, add=al)
  checkmate::assertFlag(recenter, add=al)
  
  checkmate::reportAssertions(al)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main plotting code
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  
  # make sure relevant data are available:
  x = add_cell_types_CHESS_S3(x)
  
  if (alpha_burden) {
    x = add_cell_mutation_burden_CHESS_S3(x)
  }
  
  if ("sim" %in% names(x)) {
    if (add_transformation_locations) {
      x$transformation_pos = 
        x$sim$Transformations %>% 
        dplyr::mutate(x = x + 1, y = y + 1, z = z + 1) %>%
        dplyr::mutate(label=paste0(from+1, "->", to+1))
    }
  }
      
  # center simulation
  if (recenter) {
    x = add_center_position(x)
    x = center_sim_coords(x, round(as.numeric(x$estimated_center)-c(x$params$x, x$params$y, x$params$z)))
  }
  
  # rotate simulation
  x = rotate_sim_coords(x, rotate_by, mirror)
  
  # modify coords
  type_data = x$coords
  type_data[,c("x","y","z")] = type_data[,c("x","y","z")] * mm_per_dot
  type_data$burden[type_data$type == 0] = NA
  type_data$type[type_data$type == 0] = NA
  
  if (alpha_burden) {
    # use quantiles as cut for mutation burden:
    breaks = quantile(type_data$burden, seq(0, 1, by=0.25), na.rm=1)
    type_data$burden = cut(type_data$burden, breaks, include.lowest=1, ordered=1)
  }
  
  # plot center
  if (recenter) {
    wh = x$coords$type > 0
    center = as.numeric(round(apply(x$coords[wh,c("x","y","z")], 2, mean)) * mm_per_dot)
  } else {
    if ("estimated_center" %in% names(x)) {
      center = as.numeric(x$estimated_center) * mm_per_dot
    } else {
      center = ceiling(as.numeric(x$params[c("x","y","z")]) / 2 * mm_per_dot)
    }
  }
  
  
  # construct the plot:
  space_plot = 
    ggplot(type_data, aes(x=x, y=y, fill=factor(type))) +
    cowplot::theme_cowplot()
  

  if (alpha_burden) {

    if (requireNamespace("ggrastr", quietly = TRUE)) {
      space_plot = space_plot + 
        ggrastr::geom_tile_rast(aes(alpha=burden))
    } else {
      space_plot = space_plot + 
        geom_tile(aes(alpha=burden))
    }

    space_plot = 
      space_plot + 
      labs(alpha="Mutation burden") + 
      suppressWarnings(scale_alpha_discrete(range=c(0.5, 1), na.translate=FALSE)) + # yes, using alpha for a discrete variable is not advised ...
      guides(alpha=FALSE)
    
  } else {
    
    if (requireNamespace("ggrastr", quietly = TRUE)) {
      space_plot = space_plot + 
        ggrastr::geom_tile_rast(alpha=1)
    } else {
      space_plot = space_plot +
        geom_tile(aes(alpha=1))
    }
  }
  
  space_plot = 
    space_plot +
    labs(fill="Clone") +
    geom_vline(xintercept=center[1], linetype=3, alpha=0.4) +
    geom_hline(yintercept=center[2], linetype=3, alpha=0.4) + 
    xlab("mm") + ylab("mm") + 
    theme(panel.background = element_rect(fill="transparent",colour = NA)) + 
    theme(plot.background = element_rect(fill="transparent",colour = NA)) + 
    theme(legend.position = "bottom") 
  
  
  if (!all(na.omit(type_data$type) %in% names(clone_colors))) {
    
    # fall-back option
    space_plot = 
      space_plot + 
      scale_fill_brewer(
        palette = "Set1",
        direction = -1,
        na.value = "white",
        na.translate = FALSE
      )
    
  } else {
    
    space_plot = 
      space_plot + 
      scale_fill_manual(
        breaks = names(clone_colors),
        values = clone_colors,
        na.value = "white",
        na.translate = FALSE
      )
    
  }


  # if samples in S3 object add there locations:
  if ("sample_regions" %in% names(x)) {
    for (i in seq_along(x$sample_regions$region_center)) {
      
      cent = x$sample_regions$region_center[[i]][1:2] * mm_per_dot
      diam = x$sample_regions$region_diameter[1:2] * mm_per_dot
      
      d_poly = 
        data.frame(
          x=cent[1] + diam[1] * c(-1,-1,1,1,-1), #- center[1],
          y=cent[2] + diam[2] * c(-1,1,1,-1,-1) #- center[2]
        )
      
      space_plot = space_plot + 
        geom_path(data=d_poly, aes(x=x, y=y), inherit.aes =FALSE, color="black", alpha=0.75)
    }
  }
  
  if ("samples" %in% names(x) & length(x$samples)) {
    
    polygon_paths = 
      lapply(x$samples, sample_coords_to_polygon) %>% 
      lapply("*", mm_per_dot) %>%  # scale values
      reshape::melt(measure.vars=c()) # L1 is the label
    
    sample_labels = 
      polygon_paths %>% 
      dplyr::select(L1, x=xlab, y=ylab) %>% 
      unique()
    
    if (!is.null(labeller_function)) {
      sample_labels$L1 = labeller_function(as.character(sample_labels$L1))
    }
    
    if (add_polygons) {
      space_plot = space_plot + geom_path(aes(x=x,y=y, group=L1), data=polygon_paths, inherit.aes=0)
    }
    
    space_plot = 
      space_plot + 
      ggrepel::geom_text_repel(aes(label=L1, x=x+0.25*mm_per_dot, y=y+0.25*mm_per_dot), 
                               data=sample_labels, inherit.aes=0, max.iter = 10000,
                               min.segment.length=0, force=1, size=2.5)

    
  } # end adding sample locations.
  
  if (add_transformation_locations) {
    space_plot = 
      space_plot + 
      geom_point(
        data = x$transformation_pos,
        aes(x = x * mm_per_dot, y = y * mm_per_dot),
        shape = 4,
        inherit.aes = FALSE
      )
    #geom_text_repel(data=trans_locs, aes(x=x, y=y, label=label), inherit.aes=FALSE)
  }
  
  reset_limits = recenter
  if (!is.null(x$params$max_popsize)) 
    reset_limits = reset_limits | x$params$max_popsize > 0
    
  if (reset_limits) {
    xlim = range(space_plot$data$x[which(space_plot$data$type>0)], na.rm = TRUE)
    ylim = range(space_plot$data$y[which(space_plot$data$type>0)], na.rm = TRUE)
    space_plot = space_plot + xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])
  }
  
  return(space_plot)
}


#' Wrapper function arranging multiple plots for a CHESS S3 simulation.
#'
#' @param x  A CHESS_S3 object.
#' @param target_tree (optional) tree object of class phylo to plot along the data (default: NULL).
#' @param sample_tree (optional) flag indicating if simulated sequencing data should be used for the tree reconstruction (default: false).
#' @param global_vaf_plot (optional) flag indicating if a global VAF plot should be added (default: true).
#' @param run_mobster (optional) flag indicating if mobster should be run on the global vaf data (default: false).
#' @param plot_tree_function (optional) alternative function for the tree plotting (default: NULL).
#' @param ... 
#'
#' @return A ggplot object with multiple plots arranged in a grid.
#' @export
#'
#' @examples plot.CHESS_S3(add_all_samples(new_simulation(20, seed=1), "random_cells")) # a small 20x20x1 simulation, with a tree
#' @examples plot.CHESS_S3(add_all_samples(new_simulation(20, seed=1), "random_cells"), sample_tree=TRUE) # a small 20x20x1 simulation, with a tree
#' @examples plot.CHESS_S3(add_all_samples(new_simulation(20, seed=1), "random_cells"), target_tree=ape::rtree(5)) # a small 20x20x1 simulation, with a tree
#' @examples plot.CHESS_S3(new_simulation(20, seed=1)) # a small 20x20x1 simulation
#' @examples plot.CHESS_S3(new_simulation(20, seed=1), global_vaf_plot=FALSE) # a small 20x20x1 simulation, don't plot the VAF of mutations
#' @examples plot.CHESS_S3(new_simulation(20, seed=1), run_mobster=TRUE) # a small 20x20x1 simulation, don't plot the VAF of mutations
plot.CHESS_S3 = function(x, target_tree=NULL, sample_tree=FALSE, global_vaf_plot=TRUE, run_mobster=FALSE, add_history_plot=TRUE, plot_tree_function=CHESS::plot_tree, modify_tree_function=NULL, bulk_depth = 100, min_vaf_bulk = 0.05, bins_vaf=100, n_clonal=NULL, label_drivers=TRUE, n_muts_mobster=20000, ...) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assertClass(target_tree, "phylo", null.ok = TRUE, add=al)
  checkmate::assertTRUE("phylo" %in% class(sample_tree) | (is.logical(sample_tree) & length(sample_tree) == 1), add=al)
  checkmate::assertFlag(global_vaf_plot, add=al)
  checkmate::assertFlag(run_mobster, add=al)
  checkmate::assertNumeric(n_clonal, max.len = 1, finite = TRUE, lower = 0, null.ok = TRUE, any.missing = FALSE, add=al)

  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main plotting code
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  if (!is.null(n_clonal)) {
    x$sim$ClonalMutations = n_clonal
  } else {
    x$sim$ClonalMutations = x$params$n_clonal
  }
  
  # space & history plot
  sim_plot = plot_space_CHESS_S3(x, ...) # space plot
  
  if (add_history_plot) {
    history_plot = plot_history_CHESS_S3(x) #+ guides(color=FALSE)
  } else {
    history_plot = NULL
  }
  
  
  # tree plot
  if ("phylo" %in% class(sample_tree)) {
    tree = list(tree=sample_tree)
  } else {
    tree = get_tree(x, sample_tree)
  }

  if (!is.null(tree$tree)) {
    
    # create tree plot:
    tree_umod = remove_root_tip(tree$tree, "GL")
    if (!is.null(modify_tree_function)) { tree = modify_tree_function(tree) }
    tree_mod = remove_root_tip(tree$tree, "GL")
    tree_plot = plot_tree_function(tree_mod, HI=tree$HI) # tree 
    tree_plot$labels$title = paste(tree_plot$labels$title, "(Simulation)")
    
    # optional label driver edge:
    if (label_drivers) {
      tree_plot = add_driver_labels_to_tree_plot(tree_plot, tree_umod, x) 
    }
    
  } else {
    tree_plot = NULL
  }
  
  if (!is.null(target_tree)) {
    target_tree = plot_tree_function(target_tree)
    target_tree$labels$title = paste0(target_tree$labels$title, " (Target)")
  }
  
  # global bulk sample
  if (global_vaf_plot | run_mobster) {
    mdata = get_global_vaf_sample(x, depth=bulk_depth, min_vaf=min_vaf_bulk)
  }
  
  # vaf plot
  if (global_vaf_plot) {
    global_bulk_plot = 
      plot_vaf_sample(mdata, n_bins=bins_vaf) + 
      xlab(paste0("VAF (", bulk_depth, "x)"))
  } else {
    global_bulk_plot = NULL
  }
  
  
  # mobster plot
  if (run_mobster) {
    cat("Running mobster:\n")
    
    # subset mutation data
    wh = sample(NROW(mdata), min(c(n_muts_mobster, NROW(mdata))))
    mdata_sampled = mdata[wh,]

    #
    mobster_fit =
      mdata_sampled %>%
      mutate(VAF = vaf) %>%
      mobster::mobster_fit(
        parallel = 0,
        epsilon = 1e-6,
        samples = 2,
        maxIter = 200,
        tail = TRUE
      )
    
    
    # drop tail from data
    set.seed(123)
    k_tail = mobster_fit$best$z_nk[,"Tail"]
    k_tail[k_tail < 0.01] = 0
    drop = as.logical(rbinom(seq_along(k_tail), size=1, prob=k_tail))
    d_no_tail = mobster_fit$best$data[!drop,]
    
   
    # fit mixture of binomials
    dbmix = d_no_tail[,c("alt","depth")]
    bmix_fit_b = BMix::bmixfit(dbmix, K.Binomials=1:4, samples=ifelse(NROW(d_no_tail) < 1000, 5, 2))
    bmix_params = BMix::Parameters(bmix_fit_b)
    mu_binom = bmix_params$mean %>% magrittr::set_names(bmix_params$cluster)
    
    # calculate log-lik for each variant:
    log_sum_exp = function(x) { x = x - max(x); x = exp(x); x / sum(x)}
    
    d_beta = do.call(cbind, lapply(seq_len(mobster_fit$best$K), {
        function(x) mobster::ddbpmm(x=mobster_fit$best, mdata$vaf, components=x)
      })) %>% magrittr::set_colnames(names(mobster_fit$best$N.k)) %>% 
      apply(1, log_sum_exp) %>% t()
    
    if (mobster_fit$best$K == 1) {
      d_beta = t(d_beta)
      colnames(d_beta) = names(mobster_fit$best$N.k)
    }
    
    
    d_binom = do.call(cbind, lapply(mu_binom, function(m) {
      dbinom(mdata$alt, mdata$depth, m, log = TRUE)
    })) %>% apply(1, log_sum_exp) %>% t()
    
    if (length(mu_binom) == 1) {
      d_binom = t(d_binom)
      colnames(d_binom) = names(mu_binom)
    }
    
    cluster_assignment = colnames(d_binom)[apply(d_binom, 1, which.max)]
    wh_tail = which(rbinom(length(d_beta[,"Tail"]), size=1, prob=d_beta[,"Tail"]) == 1)
    cluster_assignment[wh_tail] = "Tail"
    mdata$cluster_bmix = cluster_assignment
    
    # modify label:
    mobster_plot = 
      mdata %>% 
      mutate(clone=gsub("Bin ", "", cluster_bmix)) %>% 
      plot_vaf_sample(n_bins=bins_vaf, show_cluster_means=TRUE) + 
      xlab(paste0("VAF (", bulk_depth, "x)")) + 
      labs(fill="Cluster", color="Cluster") +
      scale_fill_viridis_d(end  = 0.8) + 
      scale_color_viridis_d(end = 0.8)
    
  } else {
    mobster_plot = NULL
  }
  
  # bind plots
  plot_list = 
    list(target_tree,
         tree_plot,
         global_bulk_plot,
         mobster_plot,
         sim_plot, 
         history_plot)
  
  plot_list = plot_list[!sapply(plot_list, is.null)]
  
  caption = 
    cowplot::ggdraw() + 
    cowplot::draw_label(parse(text=format_parameter_legend(x)), size = 10)
  
  plot_final = 
    cowplot::plot_grid(
      plotlist = plot_list,
      nrow = 1,
      align = "h",
      axis = "lbtr"
    )
  
  plot_final = 
    cowplot::plot_grid(
      plotlist = list(plot_final, caption),
      nrow = 2,
      align = "v",
      rel_heights = c(0.9, 0.1)
    )
  
  return(plot_final)
}

format_parameter_legend = function(x) {
  with(x$params, {
    
    format_param = function(x) paste0(ifelse(length(x)>1, "(", ""), paste0(signif(x, 3), collapse = "~~"), ifelse(length(x)>1, ")", ""))
    
    param_strings = apply(param.matrix, 1, format_param)
    
    symbols = c("birthrates"="lambda", "deathrates"="mu", "aggressions"=NA, 
                "push_power"="d[push]", "mutation_rates"="m", 
                "clone_start_times"="t[s]", "kill_regrow_times"=NA, 
                "father"="i[a]")
    
    param_vector = "paste(\"Parameters: \","
    for (i in seq_along(param_strings)) {
      c_symbol = symbols[names(param_strings)][i]
      if (is.na(c_symbol)) next()
      c_el = paste0(c_symbol, " == ", param_strings[i])
      param_vector = paste0(param_vector, ifelse(i>1, ",\",\",~", ""), c_el)
    }
    param_vector = paste0(param_vector, ")")
    
    return(param_vector)
  })
}

#' Function to plot a animation of a CHESS S3 simulation
#'
#' @param x CHESS S3 object. Must include snapshots (see 'snapshot_times' argument of new_simulation function (?new_simulation))
#' @param mm_per_dot (optional) value that scales the grid points to a pseudo mm value (default: 0.1).
#' @param background_color (optional) color by which to color the background of the simulated space (default: NULL, no color).
#' @param color_by  (optional) column in snapshot data by which to color cells (default: type, the cell type).
#'
#' @return A 'gganim' object containg the animation.
#' @export
#' 
#' @examples anim = animate.CHESS_S3(new_simulation(20, snapshot_times = 1:100))
animate.CHESS_S3 = function(x, mm_per_dot=0.1, background_color=NULL, color_by="type", alpha_burden=TRUE) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assert_true("snapshots" %in% names(x), add=al)
  checkmate::assertNumeric(mm_per_dot, any.missing = FALSE, null.ok = FALSE, finite = TRUE, len=1, lower=0, add=al)
  checkmate::assert_true(mm_per_dot > 0.0, add=al)
  checkmate::assert_character(background_color, len=1, any.missing = FALSE, null.ok = TRUE, add=al)
  checkmate::assert_true(color_by %in% colnames(x$snapshots$data[[1]]), add=al)
  
  checkmate::reportAssertions(al)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main plotting code
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  
  # get snapshot data from object
  snapshot_data = x$snapshots$data
  for (i in seq_along(snapshot_data)) {
    snapshot_data[[i]]$t = x$snapshots$time_exp[i]
    snapshot_data[[i]]$burden = snapshot_data[[i]]$burden / max(c(snapshot_data[[i]]$burden,1))
  }
  snapshot_data = do.call(what=rbind, snapshot_data)
  
  if (!is.null(background_color)) {
    space_polygon = with(x$params, data.frame(x=c(0, x, x, 0,0), y=c(0,0,y,y,0)))
  }
    
  animation = 
    snapshot_data %>% 
    dplyr::mutate(t=log2(t)) %>% 
    ggplot(aes(x=x*mm_per_dot, y=y*mm_per_dot))
  
  if (!is.null(background_color)) {
    animation = animation + 
      geom_polygon(data=space_polygon, aes(group=1), fill=background_color)
  }
  
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    tile_function = ggrastr::geom_tile_rast
  } else {
    tile_function = geom_tile
  }

  if (alpha_burden) {
    animation = animation + 
      tile_function(aes(fill=factor(get(color_by)), group=id, alpha=burden)) + 
      labs(alpha="Mutation burden") + 
      guides(alpha=FALSE)
  } else {
    animation = animation + 
      geom_tile(aes(alpha=burden))
  }

  
  animation = animation + 
    scale_fill_brewer(palette="Set1", direction = -1, na.value="white", na.translate=FALSE) +
    scale_alpha_continuous(range=c(0.5, 1)) + 
    guides(alpha=FALSE) +
    transition_time(t) +
    xlab("mm") + ylab("mm") + 
    theme(panel.background = element_rect(fill="transparent",colour = NA)) + 
    theme(plot.background = element_rect(fill="transparent",colour = NA)) + 
    theme(legend.position = "bottom") + 
    labs(fill="Clone") + 
    xlim(0, max(c(x$params$x, x$params$y)) * mm_per_dot) + 
    ylim(0, max(c(x$params$x, x$params$y)) * mm_per_dot) +
    labs(title = "t = {signif(frame_time, 4)} population doublings")
  
  return(animation)
}


#' Internal function to rotate grid coordinates in a data.frame
#'
#' @param d A data.frame
#' @param angle Angle in radians by which to rotate the data.
#' @param center Center of the coordinate system around which to rotate.
#' @param mirror (Optional) flag indicating if data should be mirrored around the x-axis.
#'
#' @return Data.frame with the rotated coordinates
#' @keywords internal
rotate_data = function(d, angle, center, mirror=FALSE) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  
  checkmate::assertDataFrame(d, add=al)
  checkmate::assertTRUE(all(c("x","y") %in% colnames(d)), add=al)
  checkmate::assertNumeric(angle, len=1, finite = TRUE, any.missing = FALSE, add=al)
  checkmate::assertNumeric(center, min.len=2, max.len = 3, finite = TRUE, any.missing = FALSE, null.ok = FALSE, add=al)
  checkmate::assertFlag(mirror, add=al)
  
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main function body
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  cos_a = cos(angle)
  sin_a = sin(angle)
  
  x_o = (d$x - center[1])
  y_o = (d$y - center[2])
  
  if (mirror) { 
    x_o = x_o * -1 
  }
  
  x_n = round(x_o  * cos_a - y_o * sin_a) 
  y_n = round(x_o  * sin_a + y_o * cos_a) 
  
  d$x = x_n + center[1] 
  d$y = y_n + center[2]
  
  return(d)
}


#' Internal function to offset all grid coordinates in a data.frame
#'
#' @param d A data.frame
#' @param angle Angle in radians by which to rotate the data.
#' @param center Center of the coordinate system around which to rotate.
#' @param mirror (Optional) flag indicating if data should be mirrored around the x-axis.
#'
#' @return Data.frame with the rotated coordinates
#' @keywords internal
center_data = function(d, pos) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  
  checkmate::assertDataFrame(d, add=al)
  checkmate::assertNumeric(pos, len=3, finite = TRUE, any.missing = FALSE, add=al)
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main function body
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  d$x = d$x - pos[1]
  d$y = d$y - pos[2]
  d$z = d$z - pos[3]
  
  return(d)
}


#' Internal function to rotate all coordinat data in a CHESS S3 object
#'
#' @param x CHESS S3 object.
#' @param phi (Optional) angle in radians by which to rotate the data (default: rotation angle of sample locations).
#' @param mirror (Optional) flag indicating if data should be mirrored around the x-axis (default: false).
#'
#' @return CHESS S3 object with modified coordinates.
#' @keywords internal
rotate_sim_coords = function(x, phi=NULL, mirror=FALSE) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assertNumeric(phi, any.missing = FALSE, null.ok = TRUE, finite = TRUE, len=1, add=al)
  checkmate::assertFlag(mirror > 0.0, add=al)
  checkmate::reportAssertions(al)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Rotate simulation elements
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  if (is.null(phi)) { # rotate by offset of sampling if available
    stopifnot("phi" %in% names(x)) 
    phi = -x$phi
  }
  
  if ((phi %% (2*pi) == 0) & !mirror) return(x)
  
  center = c(x$params$x, x$params$y) / 2
  
  # rotate main location data
  if ("coords" %in% names(x)) {
    x$coords = 
      x$coords %>% rotate_data(phi, center, mirror)
  }
  
  if ("transformation_pos" %in% names(x)) {
    x$transformation_pos = 
      x$transformation_pos %>% rotate_data(phi, center, mirror)
  }
  
  # rotate sample location
  for (s in names(x$samples)) {
    stopifnot(isTRUE(all.equal(x$samples[[s]]$sample_coords, x$samples[[s]]$sample_coord)))
    coords = rbind(x$samples[[s]]$sample_coords[1:3], x$samples[[s]]$sample_coords[4:6])
    coords = data.frame(coords)
    colnames(coords) = c("x","y","z")
    coords = rotate_data(coords, phi, center, mirror)
    x$samples[[s]]$sample_coord = x$samples[[s]]$sample_coords = as.numeric(c(coords[1,], coords[2,]))
  }
  
  return(x)
}


#' Internal function to rotate all coordinat data in a CHESS S3 object
#'
#' @param x CHESS S3 object.
#' @param phi (Optional) angle in radians by which to rotate the data (default: rotation angle of sample locations).
#' @param mirror (Optional) flag indicating if data should be mirrored around the x-axis (default: false).
#'
#' @return CHESS S3 object with modified coordinates.
#' @keywords internal
center_sim_coords = function(x, pos=c(0,0,0)) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assertNumeric(pos, any.missing = FALSE, null.ok = TRUE, finite = TRUE, len=3, add=al)
  checkmate::reportAssertions(al)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Shift simulation elements
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # rotate main location data
  if ("coords" %in% names(x)) {
    x$coords = x$coords %>% center_data(pos)
  }
  
  if ("transformation_pos" %in% names(x)) {
    x$transformation_pos = x$transformation_pos %>% center_data(pos)
  }
  
  
  if ("estimated_center" %in% names(x)) {
    x$estimated_center = x$estimated_center - pos
  }
  
  
  # rotate sample location
  for (s in names(x$samples)) {
    stopifnot(isTRUE(all.equal(x$samples[[s]]$sample_coords, x$samples[[s]]$sample_coord)))
    coords = rbind(x$samples[[s]]$sample_coords[1:3], x$samples[[s]]$sample_coords[4:6])
    coords = data.frame(coords)
    colnames(coords) = c("x","y","z")
    coords = center_data(coords, pos)
    x$samples[[s]]$sample_coord = x$samples[[s]]$sample_coords = as.numeric(c(coords[1,], coords[2,]))
  }
  
  return(x)
}


#' Internal function for plotting of VAF data
#'
#' @param data data.frame with columns 'vaf' and 'clone'.
#'
#' @return ggplot object.
#' @keywords internal
plot_vaf_sample = function(data, n_bins=100, show_cluster_means=FALSE) {
  
  clone_colors = c(
    "CC" = "gray25",
    "0" = "gray25",
    "1" = "#377EB7",
    "2" = "#E3211C",
    "3" = "#4DAE4A",
    "4" = "#974EA2",
    "5" = "#F680BE",
    "6" = "#FE7F00",
    "HH" = "gray50"
  )
  
  checkmate::assertDataFrame(data, col.names = "named")
  checkmate::assert_subset(c("clone","vaf"), colnames(data))
  
  # Relabel clusters:
  data = 
    data %>% 
    dplyr::mutate(clone=ifelse(clone == 0, "HH", clone)) %>% 
    dplyr::mutate(clone=ifelse(clone == -1, 0, clone))
    
  # construct title
  data_cl = data %>% select(clone, id) %>% as.data.frame() %>% unique()
  frc = table(data_cl$clone) / NROW(data_cl)
  prefix = ifelse(!is.na(suppressWarnings(as.numeric(names(frc)))), "C", "")
  clust_lab = paste0(prefix, names(frc), " = ", signif(frc*100, 2), "%")
  subtitle = paste0("N = ", NROW(data_cl), "; ", paste0(clust_lab, collapse = ", "))
  
  # calculate per cluster means (not for tail):
  if (show_cluster_means) {
    cluster_centers = 
      tapply(data$vaf, data$clone, mean) %>% 
      reshape2::melt() %>% 
      magrittr::set_colnames(c("clone","mean")) %>% 
      dplyr::filter(clone != "Tail")
  }

  plot = 
    data %>%
    dplyr::mutate(clone=factor(clone, unique(c("0","HH","CC", sort((clone)))))) %>%
    ggplot(aes(x=vaf, fill=clone, color=clone)) + 
    cowplot::theme_cowplot() + 
    geom_histogram(bins=n_bins, alpha=0.9) + 
    xlim(0, 1) +
    xlab("VAF") +
    ylab("Counts") + 
    labs(fill="Clone", color="Clone") + 
    theme(plot.title=element_text(size=9, color="gray20", face = "plain")) + 
    theme(plot.subtitle=element_text(size=8, color="gray20")) + 
    theme(legend.position = "bottom") + 
    ggtitle(subtitle)
  
  
  # set colors
  if (!all(na.omit(data$clone) %in% names(clone_colors))) {
    
    # fall-back option
    plot = 
      plot + 
      scale_fill_brewer(
        palette = "Set1",
        direction = -1,
        na.value = "white",
        na.translate = FALSE
      ) + 
      scale_color_brewer(
        palette = "Set1",
        direction = -1,
        na.value = "white",
        na.translate = FALSE
      ) 
    
  } else {
    
    plot = 
      plot + 
      scale_fill_manual(
        breaks = names(clone_colors),
        values = clone_colors,
        na.value = "white",
        na.translate = FALSE
      ) + 
      scale_color_manual(
        breaks = names(clone_colors),
        values = clone_colors,
        na.value = "white",
        na.translate = FALSE
      )
    
  }
    
  
  if (show_cluster_means) {
    plot = plot +
      geom_vline(
        data = cluster_centers,
        aes(xintercept = mean, color = as.character(clone)),
        linetype = 2, alpha = 1
      )
  }
  
  return(plot)
  
}


#' A internal function used for the plotting of trees
#'
#' @param tree object of class phylo to plot.
#' @param HI (optional) HI of the tree (unused).
#'
#' @return ggplot object
#' @export
#'
#' @examples plot_tree(ape::rtree(5))
plot_tree = function(tree, HI=NULL) {
  
  checkmate::assertClass(tree, "phylo", null.ok = FALSE)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Construct the tree
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  plt = 
    plt = tree %>% 
    ggtree::ggtree(layout="rectangular", color="gray10") +
    geom_tiplab(size=2.5, hjust=-0.25) + 
    theme(title=element_text(size=8, color="gray20")) + 
    scale_color_brewer(palette="Set1", na.value="gray30", drop=FALSE) + 
    geom_tippoint() + 
    scale_shape_manual(breaks=c("G","B","L","A"), values=c(G=16, B=15, L=18, A=8)) + 
    geom_treescale(y=-1, x=0, fontsize=2) +
    labs(color="") + 
    guides(shape=FALSE) + 
    guides(color = guide_legend(nrow=1)) +#, override.aes=aes(label="A"))) + 
    theme(legend.position="bottom") +
    theme(panel.background = element_rect(fill="transparent", colour = NA)) + 
    theme(plot.background = element_rect(fill="transparent", colour = NA))
  

  
  plt = plt + xlim(0, 1.6 * max(plt$data$x)) 
  
  return(plt)
}



#' Gets tree object from a CHESS S3 simulation object
#'
#' @param x The simulation (a CHESS S3 object).
#' @param sample_tree (optional) logical flag indicating if the underlying tree should be returned (default: true).
#' @param add_lowdepth (optional) logical flag indicating if the low depth samples should be included (default: true).
#' @param lp_depth_cutoff (optional) numeric value indicating depth below which a samples is considered to be a low depth sample  (default: 10).

#' @return A phylogenetic tree object (class phylo).
#' @export
#'
#' @examples sim = add_all_samples(new_simulation(20, seed=1), "random_cells") # a small 20x20x1 simulation
#' @examples get_tree(sim, sample_tree=TRUE) # the latent tree
#' @examples get_tree(sim, sample_tree=FALSE) # the MPR tree from the sequencing data
get_tree = function(x, sample_tree=TRUE, add_lowdepth=TRUE, lp_depth_cutoff=10, tree_manipulator=NULL, ...) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assertFlag(sample_tree, add=al)
  checkmate::assertFlag(add_lowdepth, add=al)
  checkmate::assertNumeric(lp_depth_cutoff, any.missing = FALSE, len = 1, add=al)
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main function code
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  if (is.null(x$samples)) {
    
    # 1/3) no tree
    tree = NULL
    HI = NA
    other_data = NULL
    
  } else {
    
    if (sample_tree) {
      # 2/3) direct sampling of underlying tree
      
      setup = x$samples
      
      if (!"center" %in% names(setup[[1]])) {
        setup = x$sample_list
        names(setup) = names(x$samples)
      }
      
      if (!"center" %in% names(setup[[1]])) {
        stop()
      }
      
      for (i in seq_along(setup)) {
        if ("sample_coord" %in% names(setup[[i]])) {
          setup[[i]]$sample_coord = setup[[i]]$sample_coord + 1 # shift from zero to 1 indexed
        }
      }
      
      tree =
        get_equivalent_single_cell_tree(
          x$sim,
          setup,
          depth_cutoff = lp_depth_cutoff,
          include_lowdepth = add_lowdepth,
          ...
        )
      
      
      
      # scale generation tree
      if ("obs_params" %in% names(x)) {
        set.seed(x$obs_params$sample_seed)
        if (x$obs_params$mutation_rate_scale_factor != 1) {
          
          has_singles = ape::has.singles(tree)
          if (has_singles){
            tree_o = tree
            tree = ape::collapse.singles(tree)
            if (!is.null(tree_manipulator)) tree = tree_manipulator(tree)
          }
          
          if ("dispersion" %in% rownames(x$params$param.matrix)) {
            # sample edge length from beta-binomial distribution
            disp = x$params$param.matrix["dispersion",1]
            mu = x$obs_params$mutation_rate_scale_factor
            n_gens = tree$edge.length
            tree$edge.length = rnbinom(length(n_gens), mu=mu * n_gens, size=1/disp)
          } else {
            if (!is.na(x$obs_params$mutation_rate_scale_factor)) {
              # sample edge length from poisson (sum of poisson is poisson)
              l_edge = tree$edge.length * x$obs_params$mutation_rate_scale_factor
              tree$edge.length = rpois(length(l_edge), l_edge)
            }
          }
          
          if (has_singles) {
            tree_cs = ape::collapse.singles(tree_o)
            tree_s = tree
            tree = tree_o
            tree$edge.length = rep(NA, length(tree$edge.length))
            
            
            desc_cs = Descendants(tree_cs, tree_cs$edge[,2], "tips")
            desc_s = Descendants(tree_s, tree_s$edge[,2], "tips")
            desc_o = Descendants(tree, tree$edge[,2], "tips")
            
            
            muts_per_gen_per_edge =
              lapply(seq_along(tree_s$edge.length), function(i) {
                as.numeric(table(sample(
                  seq_len(tree_cs$edge.length[i]),
                  tree_s$edge.length[i],
                  replace = TRUE
                )))
              })
            

            for (i in seq_along(tree$edge.length)) {
              wh_s_edge = which(sapply(desc_s, function(x) isTRUE(all.equal(x, desc_o[[i]]))))
              stopifnot(length(wh_s_edge) == 1)       
              
              m_edge = muts_per_gen_per_edge[[wh_s_edge]][1]
              if (is.na(m_edge)) m_edge = 0
              
              muts_per_gen_per_edge[[wh_s_edge]] = 
                muts_per_gen_per_edge[[wh_s_edge]][-1]
              
              tree$edge.length[i] = m_edge
            }
          }
        } 
      } else {
        if (!is.null(tree_manipulator)) tree = tree_manipulator(tree)
      }
      
      HI = NA
      other_data = NULL
      
    } else {
      # 3/3) simulated sequencing + max pars reconstruction
      
      tree_data = 
        phylo_from_samples(
          x, 
          add_lowdepth = add_lowdepth,
          depth_cutoff = lp_depth_cutoff,
          ...
        )
      
      tree = tree_data$tree
      if (!is.null(tree_manipulator)) tree = tree_manipulator(tree)
      phydata = tree_data$phy_data
      other_data = tree_data[!names(tree_data) %in% c("tree","phy_data")]
        
      if (is.null(phydata)) {
        HI = NA
      } else {
        HI = round(1 - phangorn::CI(tree, phydata), 2)
      }
    }
  }

  return(c(list(tree=tree, HI=HI), other_data))
}


add_driver_labels_to_tree_plot = function(plot, tree, sim) {
  
  # create annotation of highlighted edges:
  cell_types = sapply(sim$samples, function(s) {
    pos = s$sample_coords
    if (sum(pos[1:3]-pos[4:6]) > 0) return(NA)
    sim$sim$CellType(pos[1], pos[2], pos[3])
  })
  
  wh = cell_types > 1
  smps_per_clone = split(names(cell_types[wh]), cell_types[wh])
  
  edge_clone = sapply(smps_per_clone, function(x) {
    if (length(x) == 1) {
      wh_node = which(tree$tip.label == x)
    } else {
      wh_node = ape::getMRCA(x, phy = tree)
    }
    edge = which(tree$edge[,2] == wh_node)
    if (length(edge) == 0) return(0) else return(edge)
  })
  
  annot = 
    data.frame(tree$edge) %>% 
    magrittr::set_colnames(c("parent","node")) %>% 
    dplyr::mutate(driver=NA)
  
  if (length(edge_clone)) {
    annot$driver[edge_clone] = paste0("M", as.numeric(names(edge_clone))-1)
  }
  
  # add driver labels:
  plot$data = 
    merge(plot$data, annot, by=c("parent","node"), all.x = TRUE) %>% 
    as_tibble()
  
  plot =
    plot +
    ggrepel::geom_text_repel(
      aes(x = branch, label = driver),
      min.segment.length = 0,
      nudge_x = -max(ape::node.depth.edgelength(tree)) * 0.05,
      nudge_y = 0.2, 
      size = 2.5
    )
  
}


remove_root_tip = function(tree, out) {
  
  #  data = tree %>% as.treedata %>% as_tibble
  data = tidytree::as_tibble(tree)
  
  whO = which(data$label == out)
  whEP = which(data$parent == data$parent[whO] & data$node != data$node[whO] & data$parent != data$node)
  
  data$branch.length[whEP] = data$branch.length[whEP] + data$branch.length[whO]
  
  data = data[-whO,]
  
  data$parent = match(data$parent, data$node)
  data$node = match(data$node, data$node)
  
  tidytree::as.phylo(data)
}
