clonal_burden_from_tree =
  function(x) {
    wh_edge = !(x$edge[,1] %in% x$edge[,2])
    sum(x$edge.length[wh_edge])
  }


get_generator_particles = function(prior_mu, prior_push, prior_deathrate, prior_selection=NULL, prior_time=NULL, prior_sf_region_size=NULL, n_subclone=0, prior_dispersion=NULL) {

  force(prior_mu)
  force(prior_push)
  force(prior_deathrate)
  force(prior_selection)
  force(prior_time)
  force(prior_sf_region_size)
  force(prior_dispersion)
  force(n_subclone)

  stopifnot(n_subclone >= 0)
  stopifnot(n_subclone == 0 | (!is.null(prior_selection) & !is.null(prior_time) & n_subclone > 0) | (!is.null(prior_mu) & !is.null(prior_time) & n_subclone > 0))

  return(function(){

    # sample variable params
    seed = sample.int(.Machine$integer.max, 1)
    mu = prior_mu()
    push_power = prior_push()
    deathrate = prior_deathrate()
    birthrate = c(1)
    father = c(-1)
    clone_start_times = c(0)

    if (!is.null(prior_dispersion)) {
      dispersion = prior_dispersion()
    } else {
      dispersion = NULL
    }

    if (n_subclone) {

      # generate subclones parameters
      if (is.null(prior_selection)) {
        selection_coef = rep(1, n_subclone)
      } else {
        selection_coef = replicate(n_subclone, prior_selection())
      }
      start_time = sapply(selection_coef, prior_time)

      # order by start time
      ord = order(start_time)
      selection_coef = selection_coef[ord]
      start_time = start_time[ord]

      # append to previous data
      father = c(father, rep(0, length(selection_coef)))
      birthrate = c(birthrate, selection_coef)
      clone_start_times = c(clone_start_times, start_time)
    }

    #
    param.matrix =
      rbind(birthrates=birthrate,
            deathrates=deathrate,
            aggressions=1.0,
            push_power=push_power,
            mutation_rates=mu,
            clone_start_times=clone_start_times,
            kill_regrow_times=0,
            father=father,
            dispersion=dispersion)

    if (!is.null(prior_sf_region_size)) {
      sf = c(prior_sf_region_size(), rep(NA, n_subclone))
      param.matrix = rbind(param.matrix, region_scale_factor=sf)
    }

    #
    sim_params = list(param.matrix = param.matrix, seed = seed)

    return(sim_params)
  })
}


get_sample_generator = function(sample_annot, edge_distance, region_diameter, depth_lp, depth_wgs, coverage_data, default_purity=0.8, sample_width=c(G=1, B=5), alternative_offset_regions=NULL, find_center=FALSE, edge_distance_prior_strength=Inf, angle_prior_strength=Inf, ...) {

  force(sample_annot)
  force(edge_distance)
  force(region_diameter)
  force(alternative_offset_regions)
  force(find_center)
  force(edge_distance_prior_strength)
  force(angle_prior_strength)

  if (is.data.frame(sample_annot)) {

    # add purity and coverage
    barcodes = as.character(sample_annot$sample_barcode)
    sample_annot$purity = get_purity(barcodes)
    sample_annot$purity[is.na(sample_annot$purity)] = 0.8

    # add coverage
    if (!is.null(coverage_data)) {
      sample_annot$depth = coverage_data[barcodes,"coverage"]
    } else {
      sample_annot$depth = NA
    }

    sample_annot$depth[is.na(sample_annot$depth) & sample_annot$analyte == "L"] = depth_lp
    sample_annot$depth[is.na(sample_annot$depth) & sample_annot$analyte == "D"] = depth_wgs

    # add sample width
    sample_annot$diameter = sample_width[sample_annot$sample_type]
    sample_annot$center = NA

    #
    wh_row = !is.na(sample_annot$diameter)
    wh_cols = c("region","diameter","center","depth","purity")
    sample_annot = split(sample_annot[wh_row,wh_cols], seq_len(NROW(sample_annot)))
    names(sample_annot) = barcodes[wh_row]

    # convert to nested list
    sample_annot = lapply(sample_annot, as.list)
    sample_annot = lapply(sample_annot, function(d) {d$diameter = rep(d$diameter, 3); return(d)})

  }

  max_samples_per_region = max(table(sapply(sample_annot, "[[", "region")))

  all_glands = all(sapply(sample_annot, function(x) all(x$diameter == 1)))
  if (!all_glands) stop("Can't generate samples for bulks ...")


  # function to get random distance to edge
  if (edge_distance_prior_strength >= 1e6){
    .edge_distance = function() return(edge_distance)
  } else {
    .edge_distance = function() {
      sh1 = edge_distance_prior_strength * edge_distance
      sh2 = edge_distance_prior_strength * (1-edge_distance)

      d = rbeta(1, shape1 = sh1, shape2=sh2)
      while(d > 0.95 & edge_distance <= 0.95) {
        d = rbeta(1, shape1 = sh1, shape2=sh2)
      }
      return(d)
    }
  }

  # find tumour center
  if (find_center) {
    .get_center = function(x) as.numeric(round(add_center_position(x, 10000)$estimated_center))
  } else {
    .get_center = function(x) with(x$params, c(x,y,z) / 2)
  }

  if (!is.null(alternative_offset_regions)) {
    offsets_regions = alternative_offset_regions
  } else {
    offsets_regions = c(A=0, B=0.5, C=1, D=1.5)
  }

  # assume 'offsets_regions' is prior for drichlet if angle_prior_strength is small
  if (angle_prior_strength <= 1e6) {
    .offset_regions = function() {
      diff_os = c(diff(offsets_regions), 2-offsets_regions[length(offsets_regions)])
      names(diff_os) = names(offsets_regions)
      alpha = diff_os / 2 * angle_prior_strength
      angle_dri = MCMCpack::rdirichlet(1, alpha)
      offsets_regions_smp = (cumsum(angle_dri) - angle_dri[1]) * 2
      names(offsets_regions_smp) = names(offsets_regions)
      offsets_regions = offsets_regions_smp
    }
  } else {
    .offset_regions = function() { return(offsets_regions) }
  }

  # sampling function
  return(function(simulation, seed_samples=NULL) {

    # check volume of sampling regions
    if ("region_scale_factor" %in% rownames(simulation$params$param.matrix)) {
      sf = simulation$params$param.matrix["region_scale_factor", 1]
      region_diameter = ceiling(region_diameter * sf)
      if (prod(region_diameter) < max_samples_per_region) stop("Sample regions too small.\n")
    }

    # set seeds
    if (is.null(seed_samples)) seed_samples = sample(.Machine$integer.max, 1)
    set.seed(seed_samples)

    # random rotation of main region centers
    center = .get_center(simulation)
    phi_A = runif(1, 0, 2*pi)
    max_r = with(simulation$params, max(c(x,y,z) / 2))
    r_test = seq(0, max_r, by=0.25)

    # determine random angle of samples
    offset_regions = .offset_regions()

    # find closest edge
    region_center =
      lapply(offset_regions * pi, function(phi_o) {
        phi = phi_A + phi_o
        args = as.list(c(0, center - 1, phi))
        radius = ceiling(do.call(simulation$sim$GetDistanceToEdgeV2, args))

        r = radius * .edge_distance()
        x = floor(r * cos(phi) + center[1])
        y = floor(r * sin(phi) + center[2])
        z = 1

        c(x, y, z)
      })

    # list of all possible sample center positions
    .get_pos_vec = function(x)  seq(floor(-(x - 1) / 2), floor((x - 1) / 2))
    pos_locs = do.call(expand.grid, lapply(region_diameter, .get_pos_vec))

    # sample location center for each sample
    sample_regions = sapply(sample_annot, "[[", "region")
    samples_per_region = split(names(sample_annot), sample_regions)
    sample_centers = list()

    for (i in seq_along(samples_per_region)) {
      c_center = region_center[[names(samples_per_region)[i]]]

      n_needed = length(samples_per_region[[i]]) # number of positions needed
      pos_sampled = list() # list of positions sampled so far
      idx_locs_p = seq_len(NROW(pos_locs)) # still possible locations

      while (length(pos_sampled) < n_needed) {

        idxs_proposed = sample(idx_locs_p, size = n_needed, replace = FALSE)
        idx_locs_p = idx_locs_p[!idx_locs_p %in% idxs_proposed]

        offset = pos_locs[idxs_proposed, , drop = FALSE]
        centers_proposed = split(t(t(offset) + c_center), seq_len(NROW(offset)))

        .contains_cell = function(x) do.call(simulation$sim$CellType, as.list(x-1)) != 0
        is_accepted = sapply(centers_proposed, .contains_cell)

        pos_sampled = c(pos_sampled, centers_proposed[is_accepted])

      }

      pos_sampled = pos_sampled %>% magrittr::set_names(samples_per_region[[i]])
      sample_centers = c(sample_centers, pos_sampled)
    }
    centers = sample_centers %>% do.call(what = rbind)

    # create sampling setup
    sampling_setup = sample_annot
    for (i in seq_along(sampling_setup)) {
      seed = sample.int(.Machine$integer.max, 1)
      sampling_setup[[i]]$center = as.numeric(centers[names(sampling_setup)[i],])
      sampling_setup[[i]] = c(sampling_setup[[i]], min_vaf=0.0, seed=seed)
    }

    return(
      list(
        setup = sampling_setup,
        vars = list(phi = phi_A),
        seed = seed_samples,
        region_diameter = region_diameter,
        region_center = region_center
      )
    )

  })
}


labeller_function = function(x) {
  gsub("_L[0-9]+$", "*", gsub("_D[0-9]+$", "", gsub("^EPICC+_[[:alnum:]]+_", "", x)))
}
