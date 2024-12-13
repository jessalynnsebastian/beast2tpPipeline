#' @title Sample BEAST2 Trees
#'
#' @description
#' This function takes BEAST2 posterior trees and samples a subset of them,
#' making sure to keep tree names.
#'
#' @param trees_file A string of the path to the BEAST2 posterior trees file.
#' @param n_trees An integer of the number of trees to sample.
#' @param seed Optionally, a seed to use for tree sampling.
#' @param out_dir Optionally, a string of the directory to write the sampled trees to.
#'
#' @returns An `ape` multiPhylo object containing the tree sample.
#'
#' @export

sample_BEAST2_trees <- function(trees_file,
                                n_trees,
                                seed = NULL,
                                out_dir = NULL) {
  # Read all posterior trees
  trees <- ape::read.nexus(trees_file)

  # Sample n_trees trees
  if (!is.null(seed)) set.seed(seed)
  sample_indices <- sample(seq_along(trees), n_trees, replace = FALSE)
  sampled_trees <- trees[sample_indices]
  names(sampled_trees) <- paste0("tree_", sample_indices)

  return(sampled_trees)
}