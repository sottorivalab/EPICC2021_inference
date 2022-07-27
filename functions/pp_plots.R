plot_subclone_fracs = function(dset) {
  
  data = readRDS(dset)
  title = paste0(basename(dirname(dirname(dset))), " - ", format_model_string(basename(dirname(dset))))
  
  scs = sort(unique(unlist(data$per_sample_sc_freq[,-c(1:5)])))
  
  sc_per_sample = 
    lapply(scs, function(x) {
      scf = apply(data$per_sample_sc_freq[,-c(1:5)] == x, 2, mean)
      data.frame(sample=names(scf), frac=scf, subclone=x)
    }) %>% do.call(what=rbind)
  
  sc_per_sample = cbind(sc_per_sample, THmisc::annotation_from_barcode(sc_per_sample$sample))
  sc_per_sample$sample = gsub("EPICC_", "", sc_per_sample$sample)
  sc_per_sample$sample = factor(sc_per_sample$sample, unique(rev(sort(sc_per_sample$sample))))
  
  plot = 
    sc_per_sample %>% 
    ggplot(aes(x=sample, y=factor(subclone), fill=frac)) + 
    geom_tile() + 
    coord_flip() + 
    xlab("") + 
    ylab("Clone") + 
    facet_grid(region~., scales = "free") + 
    scale_fill_viridis_c(begin=0.2, end=0.8) + 
    geom_text(aes(label=signif(frac, 3))) + 
    theme(strip.text.y =element_text(angle=0)) + 
    labs(fill="Frac.") + 
    ggtitle(title)
  
}



plot_region_positions = function(dset) {
  
  data = readRDS(dset)
  stopifnot("sample_region_positions" %in% names(data))
  title = paste0(basename(dirname(dirname(dset))), " - ", format_model_string(basename(dirname(dset))))
  regions_to_plot = unique(data$sample_dist_data$R1)
  
  
  # priors for sampling:
  prior_rel_dist = get(".edge_distance", envir = environment(data$variables$s_gen))
  
  angle_prior_strength = get("angle_prior_strength", envir = environment(data$variables$s_gen))
  offsets_regions = get("alternative_offset_regions", envir = environment(data$variables$s_gen))
  if(is.null(offsets_regions)) offsets_regions = c(A=0, B=0.5, C=1, D=1.5)
  
  angle_prior = function() {
    diff_os = c(diff(offsets_regions), 2-offsets_regions[length(offsets_regions)])
    names(diff_os) = names(offsets_regions)
    alpha = diff_os / 2 * angle_prior_strength
    angle_dri = MCMCpack::rdirichlet(1, alpha)
    offsets_regions_smp = (cumsum(angle_dri) - angle_dri[1]) * 2
    names(offsets_regions_smp) = names(offsets_regions)
    offsets_regions_smp * pi
  }
  
  
  # sample priors
  d_prior_rel_dist = 
    data.frame(
      region="Prior", 
      rel_d_center_region = sapply(1:2000, function(x) prior_rel_dist())
    )
  
  d_prior_angles = 
    reshape::melt(do.call(rbind, lapply(1:2000, function(x) angle_prior()))) %>% 
    dplyr::mutate(type="Prior") %>% 
    dplyr::select(region=X2, angle=value, type)
  
  
  plot_rel_distance = 
    data$sample_region_positions %>% 
    dplyr::filter(region %in% regions_to_plot) %>%
    dplyr::select(region, rel_d_center_region) %>% 
    rbind(d_prior_rel_dist) %>% 
    dplyr::mutate(region = factor(region, c("Prior", LETTERS))) %>% 
    dplyr::mutate(rel_d_center_region=scales::squish(rel_d_center_region, range(0, max(d_prior_rel_dist$rel_d_center_region)))) %>% 
    ggplot(aes(x=region, y=rel_d_center_region)) + 
    ggsignif::geom_signif(comparisons =lapply(regions_to_plot, function(x) c("Prior", x)), step_increase = 0.25) +
    THmisc::pretty_violin() + 
    xlab("Region") + 
    ylab("Relative distance to center")+ 
    cowplot::theme_cowplot() + 
    scale_y_continuous(expand = expansion(mult=c(0, 0.1)), limits = c(0, NA))
  
  
  plot_distance = 
    data$sample_region_positions %>% 
    dplyr::filter(region %in% regions_to_plot) %>%
    ggplot(aes(x=region, y=d_center_region)) + 
    THmisc::pretty_violin() + 
    xlab("Region") + 
    ylab("Distance to center")+ 
    cowplot::theme_cowplot()
  
  
  plot_angle = 
    data$sample_region_positions %>% 
    dplyr::mutate(type="Post.") %>%
    dplyr::select(region, angle, type) %>% 
    dplyr::mutate(type=factor(type, c("Prior","Post."))) %>% 
    rbind(d_prior_angles) %>%
    dplyr::filter(region != "A") %>%
    dplyr::filter(region %in% regions_to_plot) %>%
    ggplot(aes(x=type, y=angle/(2*pi)*360)) + 
    ggsignif::geom_signif(comparisons = list(c("Prior", "Post."))) +
    THmisc::pretty_violin() + 
    facet_wrap(~region) + 
    xlab("Region") + 
    ylab("Relative angle") + 
    cowplot::theme_cowplot() + 
    scale_y_continuous(expand = expansion(mult=c(0, 0.1)), limits = c(0, NA))
  
  
  mean_pos_in_regions = 
    dplyr::group_by(data$sample_region_positions, region) %>% 
    dplyr::summarise(angle=median(angle), rel_d_center_region=median(rel_d_center_region))
  
  plot_example = 
    data$sample_region_positions %>% 
    dplyr::filter(region %in% regions_to_plot) %>%
    ggplot(aes(x=angle/(2*pi)*360, y=rel_d_center_region,  color=region)) + 
    geom_point(alpha=0.5, size=0.5) + 
    geom_point(data=mean_pos_in_regions, size=4, shape="x", color="gray10") + 
    xlab("Angle") + 
    ylab("Relative distance to center") +
    coord_polar(direction = -1) + 
    scale_x_continuous(breaks=c(0, 90, 180, 270), limits=c(0, 360)) +
    ylim(0, 1) + 
    scale_color_brewer(palette = "Set1") + 
    labs(color="Region") + 
    cowplot::theme_cowplot() + 
    background_grid()
  
  title_gr = ggdraw() + 
    draw_label(
      title,
      fontface = 'bold', x = 0, hjust = 0
    )
  
  
  plot_grid = 
    cowplot::plot_grid(
      plot_rel_distance,
      plot_distance,
      plot_angle,
      plot_example
    )
  
  plot_grid = 
    cowplot::plot_grid(
      title_gr,
      plot_grid, 
      rel_heights = c(0.1, 1), 
      ncol=1
    )
  
  plot_grid
}
