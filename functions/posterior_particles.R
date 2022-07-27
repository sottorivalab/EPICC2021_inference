
load_posterior_particles = function(x, last_n_states=Inf) {
  
  stopifnot(file.exists(x))
  stopifnot(is.numeric(last_n_states) & length(last_n_states) ==  1)
  
  # load data
  print(x)
  d = readRDS(x)
  
  # Keep obs below eps and from last n states:
  wh = d$observations$dist <= d$states$epsilon[which.max(d$states$i)] &
    d$observations$i %in% tail(sort(d$states$i), last_n_states)
  d$observations = d$observations[which(wh),]
  
  CHESS::drop_unobserved(d)
}


add_sampled_information = function(x, n_sample_clones, n_cores=1) {
  
  if (n_sample_clones == 0) return(x)
  
  # sample a subset of simulations:
  set.seed(123)
  n = min(c(n_sample_clones, NROW(x$observations)))
  idx = sample(NROW(x$observations), n)
  
  # iterate over simulations:
  sample_infos = pbmcapply::pbmclapply(idx, function(i) {
    
    # generate simulation
    idx_sim = unlist(x$observations[i,letters[9:12]])
    sim = CHESS::get_sim_by_idx(x, idx_sim, add_samples = TRUE)
    sim_info = x$observations[i,1:5]
    
    
    # get sample swap infos:
    target = x$variables$target_y
    rho = x$variables$rho
    dist_exp = as.numeric(x$observations[i,"dist"])
    tree_distance_mod = function(...) CHESS::tree_distance(..., return_tree=TRUE)
    assign("tree_distance", tree_distance_mod, envir=environment(rho)) # patch distance function
    
    tree = CHESS::get_tree(sim)
    tree2 = rho(x$variables$target_y, tree$tree)
    if(tree2$dist != dist_exp) return(NULL)
    stopifnot(isTRUE(dplyr::all_equal(tree$tree$edge, tree2$tree$edge)))
    stopifnot(isTRUE(dplyr::all_equal(tree$tree$edge, tree2$tree$edge)))
    sample_id_map = tree2$tree$tip.label %>% magrittr::set_names(tree$tree$tip.label)
    
    # modify sample names:
    names(sim$sample_list) = sample_id_map[names(sim$sample_list)]
    
    # get genetic distances:
    genetic_dists = cophenetic(tree2$tree)
    
    # get sample positions:
    sample_centers = do.call(rbind, lapply(sim$sample_list, "[[", "center"))
    sample_dists = as.matrix(dist(sample_centers))
    
    # get subclone types:
    if ("clone_start_times.2" %in% colnames(x$particles)) {
      subclone_type =
        sapply(sim$sample_list, function(s) {
          sim$sim$CellType(s$center[1]-1, s$center[2]-1, s$center[3]-1)
        }) %>% data.frame %>% t()
    } else {
      subclone_type = rep(1, length(sim$sample_list))
    }
    
    # get sample region center
    sample_seed = CHESS::get_sample_seed(x, idx_sim)
    sample_setup = x$variables$s_gen(sim, sample_seed)
    sim_center = as.numeric(CHESS::add_center_position(sim)$estimated_center)
    
    .get_dist = function(x, y) {as.numeric(dist(rbind(x,y)))}
    .get_angle = function(x, y){ atan2(as.numeric(x[2]-y[2]), as.numeric(x[1]-y[1]))}
    .dist_to_edge = function(pos, angle) sim$sim$GetDistanceToEdgeV2(0, pos[1]-1, pos[2]-1, pos[3]-1, angle)
    
    dist_to_center = sapply(sample_setup$region_center, .get_dist, y=sim_center)
    angle = sapply(sample_setup$region_center, .get_angle, y=sim_center)
    dist_to_edge = sapply(angle, .dist_to_edge, pos=as.numeric(sim_center))
    rel_pos = dist_to_center / dist_to_edge

    region_pos_data =
      data.frame(
        region = names(angle), 
        angle = (angle-angle[1])%%(2*pi),
        d_center_edge = dist_to_edge, 
        d_center_region = dist_to_center, 
        rel_d_center_region = rel_pos
      )
    
    
    return(list(sc_data=cbind(sim_info, subclone_type), dists=sample_dists, gdists=genetic_dists, region_pos_data=region_pos_data))
  }, mc.cores = n_cores)
  
  # drop failed sims
  wh_bad = sapply(sample_infos, is.null)
  if (any(wh_bad)) warning(paste0("Couldn't reproduce ", signif(mean(wh_bad)*100, 3), "% of simulations."))
  sample_infos = sample_infos[!wh_bad]

  smp_ids = sort(rownames(sample_infos[[1]]$dists))
  
  x$per_sample_sc_freq = 
    lapply(sample_infos, "[[", "sc_data") %>% 
    do.call(what=rbind) %>% 
    tibble::as_tibble()
  
  sample_dists = lapply(sample_infos, function(x) x$dists[smp_ids,smp_ids])
  
  x$sample_dist_data = 
    do.call(abind::abind,  list(sample_dists, along=3)) %>% 
    reshape2::melt() %>% magrittr::set_colnames(c("S1","S2","i","dist")) %>% 
    dplyr::filter(as.character(S1) < as.character(S2)) %>% 
    dplyr::mutate(R1 = substr(S1, 12, 12)) %>% 
    dplyr::mutate(R2 = substr(S2, 12, 12))
  
  
  sample_gdists = lapply(sample_infos, function(x) x$gdists[smp_ids,smp_ids])
  
  x$sample_genetic_dist_data = 
    do.call(abind::abind,  list(sample_gdists, along=3)) %>% 
    reshape2::melt() %>% magrittr::set_colnames(c("S1","S2","i","dist")) %>% 
    dplyr::filter(as.character(S1) < as.character(S2)) %>% 
    dplyr::mutate(R1 = substr(S1, 12, 12)) %>% 
    dplyr::mutate(R2 = substr(S2, 12, 12))
  
  
  x$sample_region_positions = 
    lapply(sample_infos, "[[", "region_pos_data") %>% 
    do.call(what=rbind) %>% 
    tibble::as_tibble()
  
  x
}

