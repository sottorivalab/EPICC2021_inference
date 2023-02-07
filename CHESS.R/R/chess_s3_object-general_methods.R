#' Constructor for a CHESS simulation S3 object (CHESS_S3).
#'
#' @param x Size of dimension x. Soft boundary if set to 1.
#' @param y Optional size of second dimension (default: same as x).
#' @param z Optional size of third dimension (default: 1)..
#' @param param.matrix 
#' @param seed 
#' @param verbose 
#' @param clonal_mutations 
#' @param snapshot_times 
#' @param explore_locally 
#' @param record_history 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples new_simulation(x=20, y=20) # 20x20 simulation
new_simulation = # CHESS_S3 const
  function(x, y=x, z=1, param.matrix=NULL, seed=NULL, verbose=FALSE, clonal_mutations=0, snapshot_times=NULL, explore_locally=FALSE, record_history=TRUE, return_generation_tree=FALSE, max_popsize=0, n_center_estimate=0, alt_edge_finding=FALSE, ...) {
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Check arguments
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    al = checkmate::makeAssertCollection()
    
    # dimensions
    checkmate::assertIntegerish(x, any.missing = FALSE, null.ok = FALSE, len = 1, lower=1, add=al)
    checkmate::assertIntegerish(y, any.missing = FALSE, null.ok = FALSE, len = 1, lower=1, add=al)
    checkmate::assertIntegerish(z, any.missing = FALSE, null.ok = FALSE, len = 1, lower=1, add=al)
    checkmate::assert_true(max(c(x, y, z)) >= 3)
    
    
    # test param.matrix
    checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE, null.ok = TRUE, add=al)
    checkmate::assertFlag(verbose, add=al)
    checkmate::assertIntegerish(clonal_mutations, len = 1, lower=0, add=al)
    checkmate::assertNumeric(snapshot_times, lower=1, any.missing = FALSE, min.len = 0, finite = TRUE, null.ok = TRUE, add=al)
    checkmate::assertNumeric(n_center_estimate, lower=0, any.missing = FALSE, min.len = 0, finite = TRUE, null.ok = TRUE, add=al)
    checkmate::assertFlag(explore_locally, add=al)
    checkmate::assertFlag(record_history, add=al)
    checkmate::assertFlag(return_generation_tree, add=al)
    
    checkmate::reportAssertions(al)
    
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Set default arguments (optional)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
    # default is a random seed
    if (is.null(seed)) {
      seed = sample.int(.Machine$integer.max, 1)
    }
    
    param.matrix_default =
      rbind(
        birthrates = 1,
        deathrates = 0,
        aggressions = 1.0,
        push_power = max(c(x, y, z)),
        mutation_rates = 10,
        clone_start_times = 0,
        kill_regrow_times = 0,
        father = -1
      )
    
    # as default neutral simulation
    if (is.null(param.matrix)) {
      param.matrix = param.matrix_default
    }
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Generate the simulation
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    # reorder subclones based on start time
    ord = order(param.matrix["clone_start_times",])
    param.matrix = param.matrix[,ord,drop=FALSE]

    # catch NA mutation rates:
    if (all(is.na(param.matrix["mutation_rates",]))) {
      param.matrix["mutation_rates",] = 0
      return_generation_tree = TRUE
    }
    
    
    # run simulation with passed arguments
    param.matrix_use = param.matrix[rownames(param.matrix_default),,drop=FALSE]
    
    flags = list(
      explore_locally = explore_locally,
      return_generation_tree = return_generation_tree,
      max_popsize = max_popsize,
      n_center_estimate = n_center_estimate
    )
    
    params = list(
      x = x,
      y = y,
      z = z,
      param.matrix = param.matrix_use,
      seed = seed,
      n_clonal = clonal_mutations
    )
    
    sim = do.call(new, c(CHESS_Universe_object, params))
    params$param.matrix = param.matrix
    
    # set flags and optional parameters:
    if (length(snapshot_times)) sim$SnapshotTimes = as.numeric(snapshot_times)
    if (explore_locally) sim$ExploreLocally  = TRUE 
    if (record_history) sim$RecordHistory = TRUE
    if (return_generation_tree) sim$GenerationTree = TRUE
    if (max_popsize > 0) sim$StopPopulationSize = max_popsize
    if (alt_edge_finding) sim$AltPush = TRUE
    
    # execute simulation
    finished = sim$RunSimulation(verbose)
    if (!finished) {
      warning("Simulation was interrupted. Not returning results.")
      return(invisible(NULL))
    }
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Construct results
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    # construct the S3 object:
    res = list(params=c(params, flags), sim=sim, is_finished=finished)
    class(res) = "CHESS_S3"
    
    # add snapshot data
    if (length(snapshot_times)) {
      
      snapshot_data = sim$Snapshots
      
      res$snapshots = 
        list(
          data = snapshot_data,
          time = sim$SnapshotTimes,
          time_exp = snapshot_times[seq_along(snapshot_data)]
        )
      
    }
    
    # estimate tumour center
    if (n_center_estimate > 0) {
      res = add_center_position(res, n_center_estimate)
    }
    
    
    return(res)
  }


print.CHESS_S3 = function(x, ...) {

  cat("A CHESS simulation\n")
  
  if ("is_finished" %in% names(x)) {
    if (!x$is_finished) {
      cat(crayon::bgRed("**** Simulation is unfinished! ****\n"))
    }
  }
  
  dim_string = paste(x$params[c("x","y","z")], collapse=" x ")
  cat(" - Dimensions:", dim_string, "\n")
  cat(" - Clonal mutations:", x$params$n_clonal, "\n")
  cat(" - Seed:", x$params$seed, "\n")
  cat(" - Clones:\n")
  print(data.frame(t(x$params$param.matrix)))
  cat("\n\n")
  cat("Available methods:\n")
  cat(" - print(x)\n", sep="")
  cat(" - plot(x)\n", sep="")
  if ("snapshots" %in% names(x)) {
    cat(" - animate(x)\n", sep="")
  }
}


phyDat_to_tree = function(phyDat, all_trees=FALSE) {
  
  phyDat %>%
    phangorn::pratchet(trace=FALSE, all=all_trees) %>%
    ape::unroot() %>%
    ape::root(out="GL", resolve.root=TRUE) %>%
    phangorn::acctran(phyDat)
  #remove_root_tip("GL")
}


get_tips_below = function(i, tree) {
  
  if (i %in% tree$tip.label) 
    return(i)
  
  if (!is.numeric(i))
    i = as.numeric(as.character(i))
  
  stopifnot(!is.na(i))
  
  if (i <= length(tree$tip.label)) 
    return(tree$tip.label[i])
  
  edge_below = tree$edge[tree$edge[,1]==i, 2] 
  if (length(edge_below) == 0) return(edge_below)
  
  unlist(sapply(edge_below, get_tips_below, tree=tree))
}
