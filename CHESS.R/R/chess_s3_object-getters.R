#' Get phylogenetic data from samples in a CHESS S3 object
#'
#' @param x A CHESS_S3 object. Must include samples. 
#' @param all_trees (optional) flag indicating whether all maximum parsimony trees should be returned (default: false).
#' @param add_lowdepth (optional) flag indicating whether low coverage samples should be added through a ML method (default: false).
#' @param min_confidence (optional) numeric value indicating minimal confidence for the ML method (default: 0.99).
#' @param ... (optional) arguments passed to the phy_data_from_samples or add_lowpass_sampled function (see ?phy_data_from_samples and ?add_lowpass_sampled for details).
#'
#' @return A phylogentic tree.
#' @export
#'
phylo_from_samples = function(x, all_trees=FALSE, add_lowdepth=FALSE, min_confidence=0, ...) {
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  al = checkmate::makeAssertCollection()
  checkmate::assertClass(x, "CHESS_S3", add=al)
  checkmate::assert_true("samples" %in% names(x), add=al)
  checkmate::assertFlag(all_trees, add=al)
  checkmate::assertFlag(add_lowdepth, add=al)
  checkmate::assertNumeric(min_confidence, lower = 0, upper = 1, len = 1, any.missing = FALSE, add=al)
  checkmate::reportAssertions(al)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main code
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # construct tree
  phy_data = phy_data_from_samples(x, ...)
  tree = phyDat_to_tree(phy_data, all_trees=all_trees)
  lp_adding_data = NULL
  
  # optionally, assign low pass data to tree:
  if (add_lowdepth) {
    
    wh_missing = !names(x$samples) %in% tree$tip.label
    lp_data = lapply(x$samples[wh_missing], "[[", "mutation_data")
    
    if (length(lp_data)) { # any lp samples?
      
      lp_adding_data = 
        MLLPT::add_lowpass_sampled(
          tree = tree,
          phydata = phy_data, 
          sample_data = lp_data,
          min_confidence = min_confidence,
          return_details = FALSE,
          ...
        )
      
      tree = lp_adding_data
      phy_data = NULL

    }
  } 
  
  result_list = list(tree=tree, phy_data=phy_data)
  
  return(result_list)
}



phy_data_from_samples = function(d, ccf_cutoff = 0.5, depth_cutoff=10, include_lowdepth=FALSE, ...){
  
  # identify lowpass samples
  mean_depth = sapply(d$samples, function(x) mean(x$mutation_data$depth))
  wh_use = mean_depth >= depth_cutoff
  
  # create copy excluding lowpass samples
  d_cp = d[c("params","samples")]
  d_cp$sim = 0
  d_cp$samples[!wh_use] = NULL
  
  # get vaf data
  vaf_table = vaf_table_from_samples(d_cp)
  pur = sapply(d$samples, "[[", "purity")[colnames(vaf_table)]
  ccf_table = data.frame(t(t(vaf_table) / pur * 2))
  ccf_table$GL = 0
  
  # to phydat
  mutated = ccf_table > ccf_cutoff
  storage.mode(mutated) = "numeric"
  pyData = mutated %>% t() %>% phangorn::phyDat(type="USER", levels=c(0,1))
  attr(pyData, "id") = rownames(mutated)
  
  return(pyData)
}



table_from_samples = function(d, value="vaf", min_value=NULL) {
  
  if (!"samples" %in% names(d)) { return(NULL) }
  stopifnot(value %in% names(d$samples[[1]]$mutation_data))
  
  samples = d$samples
  all_ids = unique(unlist(lapply(samples, function(x) x$mutation_data$id)))
  mdata_all = data.frame(variant = all_ids)
  
  for (i in seq_along(samples)) {
    sample_id = names(samples)[i]
    mdata_sample = samples[[i]]$mutation_data
    mdata_all[,sample_id] = 0 # default vaf 0.0
    mt = match(mdata_sample$id, mdata_all$variant)
    mdata_all[mt,sample_id] = mdata_sample[,value]
  }
  
  # drop variant column
  rownames(mdata_all) = mdata_all$variant
  mdata_all$variant = NULL
  
  # filter based on max value
  if (!is.null(min_value)) {
    max_value = apply(mdata_all, 1, max)
    mdata_all = mdata_all[max_value >= min_value,]
  }
  
  return(mdata_all)
  
}

vaf_table_from_samples = function(d, min_vaf_filter=0.05) {
  table_from_samples(d, value="vaf", min_value = min_vaf_filter)
}


get_equivalent_single_cell_tree = function(s, setup, depth_cutoff=10, include_lowdepth=TRUE) {
  
  stopifnot("Rcpp_CHESS_Universe_object" %in% class(s))
  stopifnot(is.logical(include_lowdepth) & is.numeric(depth_cutoff))
  stopifnot(length(include_lowdepth) == 1 & length(depth_cutoff) == 1)
  
  # optionally remove low depth samples
  is_low_depth = sapply(setup, function(x) x$depth < depth_cutoff)
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
  
  # get the tree
  single_cell_tree_data = s$SingleCellTree(x_cells, y_cells, z_cells)
  tree_object = ape::read.tree(text=single_cell_tree_data$tree_string)
  tree_object = ape::collapse.singles(tree_object) # collapse all single nodes
  
  
  # replace labels with the ones in the setup object
  coord_label_tree = gsub(".*_type[[:digit:]]+_", "", tree_object$tip.label)
  tree_object$tip.label = names(setup)[match(coord_label_tree, coord_label)]
  tree_object = ape::root.phylo(add.tips(tree_object, "GL", getRoot(tree_object), 0), "GL")
  
  
  # reassign lowpass samples
  if (any(is_low_depth)) {
    
    min_edge_length = floor(max(tree_object$edge.length)*0.025)
    label_dropped = names(which(is_low_depth))
    idx_dropped = which(tree_object$tip.label %in% label_dropped)
    
    for (i in idx_dropped) {
      
      c_edge = tree_object$edge[,2] == i # terminal edge 
      c_node = tree_object$edge[c_edge,1]
      
      while (all(get_tips_below(c_node, tree_object) %in% label_dropped)) {
        c_edge = tree_object$edge[,2] == c_node # terminal edge 
        c_node = tree_object$edge[c_edge,1]
      }
      
      new_label =  paste0(tree_object$tip.label[i], "added")
      stopifnot(!new_label %in% tree_object$tip.label)
      tree_object = phangorn::add.tips(tree_object, new_label, c_node, min_edge_length)
    }
    
    tree_object = ape::drop.tip(tree_object, label_dropped)
    tree_object$tip.label = gsub("added", "", tree_object$tip.label)
  }
  
  return(tree_object) 
}


#' Gets a matrix of VAF from samples in the CHESS S3 object.
#'
#' @param x CHESS S3 object.
#' @param ... (Optional) addtional arguments passed to the add_sample method (see ?add_sample).
#'
#' @return Returns data.frame with mutation information.
#' @export
#'
#' @examples get_global_vaf_sample(new_simulation(20, seed=1), depth=50) # a small 20x20x1 simulation
get_global_vaf_sample = function(x, ...) {
  checkmate::assertClass(x, "CHESS_S3")
  center = ceiling(unlist(x$params[c("x","y","z")]) / 2)
  diameter = unlist(x$params[c("x","y","z")])
  x = add_sample(x, center = center, diameter = diameter, ...)
  mdata = x$samples[[length(x$samples)]]$mutation_data
  
  # mark clonal cluster:
  if (x$sim$ClonalMutations) {
    cl_ids = unique(sapply(strsplit(mdata$id, "X"), "[[", 1))
    cl_ids_n = strtoi(cl_ids, 16)
    regex_clonal = paste0("^", cl_ids[which.min(cl_ids_n)], "X")
    mdata$clone = ifelse(grepl(regex_clonal, mdata$id), -1, mdata$clone)
  }

  return(mdata)
}
