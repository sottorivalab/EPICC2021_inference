plot_tree = function(tree, annot_data_tree=NULL, tree_data=NULL, color_by=NULL, size_edge_by=NULL, color_edge_by=NULL, driver_labels=NULL, HI=NULL, tip_label_size=2.5, lty_lp=1, alpha_lp=1, ...) {
  
  tree = MLLPT::set_lp_tiplength(tree, 0)
  
  require(magrittr)
  require(ggtree)
  
  #wh_not_GL = tree$tip.label != "GL"
  tip_labels_orig = tree$tip.label
  rownames(tree$edge) = seq_len(NROW(tree$edge))
  tree = ape::drop.tip(tree, "GL")
  annotation = annotation_from_barcode(tree$tip.label, extract=TRUE)
  added_node = grepl("Added", tree$tip.label) | annotation$analyte != "D"
  
  new_tip_label = gsub("EPICC_C[[:digit:]]*_", "", annotation$tissue_barcode)
  new_tip_label[added_node] = paste0(new_tip_label[added_node], "*")
  tree$tip.label = new_tip_label
  
  shape_lookup_table = annotation$sample_type %>% set_names(new_tip_label)
  shape_lookup_table[added_node] = "A"
  
  if (is.null(color_by)) {
    color_lookup_table = annotation$region %>% set_names(new_tip_label)
    color_lookup_table = factor(color_lookup_table, levels=LETTERS[1:8])
  } else {
    id = as.character(annotation$tissue_barcode)
    color_lookup_table = color_by[id] %>% set_names(new_tip_label)
  }
  
  if (is.null(size_edge_by)) {
    if (is.null(color_edge_by)) {
      plt = tree %>% 
        ggtree::ggtree(layout="rectangular", color="gray10", size=0.5, aes(alpha=alpha_edge, linetype=lty))
    } else {
      plt = tree %>% 
        ggtree::ggtree(layout="rectangular", aes(color=get(color_edge_by), alpha=alpha_edge, linetype=lty)) +
        labs(color=color_edge_by)
    }
  } else {
    if (is.null(color_edge_by)) {
      plt = tree %>% 
        ggtree::ggtree(layout="rectangular", color="gray10", size=0.5, aes(alpha=alpha_edge, size=get(size_edge_by), linetype=lty))
    } else {
      plt = tree %>% 
        ggtree::ggtree(layout="rectangular", aes(color=get(color_edge_by), size=get(size_edge_by), alpha=alpha_edge, linetype=lty)) + 
        guides(size=FALSE)
    }
  }
  
  plt$data$alpha_edge = 1
  plt$data$lty = 1
  #plt$data$color = "gray10"
  
  
  # set linetype of lp samples:
  get_tips = function(i) {
    wh = which(plt$data$parent == i)
    wh = wh[wh != i]
    if (length(wh) == 0) return(plt$data$label[which(plt$data$node == i)])
    unlist(sapply(wh, get_tips))
  }
  
  all_lp_desc = function(i) { all(grepl("[*]", get_tips(i)))}
  wh_lp_desc = sapply(plt$data$node, all_lp_desc)
  plt$data$alpha_edge[wh_lp_desc] = alpha_lp 
  plt$data$lty[wh_lp_desc] = lty_lp
  #plt$data$color[wh_lp_desc] = "gray10"
  
  if (!is.null(annot_data_tree)) {
    stopifnot(all(c("parent","node") %in% colnames(annot_data_tree)))
    plt$data = tibble(merge(plt$data, annot_data_tree, by=c("parent","node"), all.x = TRUE))
  }
  
  # 
  mt = match(plt$data$node, tree$edge[,2])
  for (el in names(tree)) {
    if (el %in% c("edge","Nnode","tip.label","edge.length")) next()
    stopifnot(length(tree[[el]]) != length(mt))
    plt$data[,el] = tree[[el]][as.numeric(rownames(tree$edge)[mt])]
  }
  
  
  if ((!is.null(color_edge_by) | !is.null(size_edge_by)) & is.null(annot_data_tree)) {
    plt = plt + 
      scale_color_continuous() + 
      scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.5, limits=c(0, 1)) + 
      ggnewscale::new_scale("color") 
  }
  
  
  plt = plt + 
    geom_tiplab(aes(color=color_lookup_table[label]), size=tip_label_size, hjust=-0.2, alpha=1) + 
    ggtitle(annotation$patient[1]) + 
    theme(title=element_text(size=8, color="gray20")) + 
    scale_color_brewer(palette="Set1", na.value="gray30", drop=FALSE) + 
    geom_tippoint(aes(shape=shape_lookup_table[label])) + 
    scale_shape_manual(breaks=c("G","B","L","A"), values=c(G=16, B=15, L=18, A=8)) + 
    geom_treescale(y=-1, x=0, fontsize=2) +
    labs(color="") + 
    guides(shape=FALSE) + 
    guides(color = guide_legend(nrow=1)) +#, override.aes=aes(label="A"))) + 
    theme(legend.position="bottom") +
    theme(panel.background = element_rect(fill="transparent", colour = NA)) + 
    theme(plot.background = element_rect(fill="transparent", colour = NA)) + 
    scale_alpha_identity() + 
    scale_linetype_identity()
  

  if (is.null(color_by)) { 
    plt = plt + labs(color="Region") + guides(color=FALSE)
  }
  
  if ("edge.label" %in% names(tree)) {
    stopifnot(length(tree$edge.length) == NROW(tree$edge))
    wh_label = tree$edge.label != ""
    if (sum(wh_label)) {
      edge_data = data.frame(tree$edge, edge_label=tree$edge.label)[wh_label,]
      seed = sample(1:1e5, 1)
      plt = plt %<+% edge_data + ggrepel::geom_label_repel(aes(x=branch, label=edge_label), label.size=NA, min.segment.length=0, max.iter = 10000, segment.alpha=1, nudge_x = 2, force=150, segment.size = 0.3, color="gray10", segment.color = "gray20", hjust=0, vjust=1, lineheight=.75, alpha=0.6, size=1.5, seed=seed, box.padding = 0, label.r = 0) + ggrepel::geom_label_repel(aes(x=branch, label=edge_label), label.size=NA, min.segment.length=0, max.iter = 10000, segment.alpha=1, nudge_x = 2, force=150, segment.size = 0.3, color="gray10", segment.color = "gray20", hjust=0, vjust=1, lineheight=.75, alpha=1.0, size=1.5, seed=seed, box.padding = 0, fill=NA, label.r = 0)
    }
  }
  
  if (!is.null(HI)) {
    max_x = max(plt$data$x)
    HI_color = c("TRUE"="#f9a800","FALSE"="gray10")[as.character(HI>0.1)]
    plt = plt + geom_text(data=data.frame(x=max_x, y=-1), aes(x=x*0.8, y=y), 
                          label=paste0("HI = ", HI), size=2.5, 
                          inherit.aes=FALSE, hjust=0, color=HI_color)
  }
  
  plt = plt + xlim(0, 1.6 * max(plt$data$x)) 
}

plot_tree_alt = plot_tree

