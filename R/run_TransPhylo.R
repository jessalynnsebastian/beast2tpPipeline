#' @title Run TransPhylo on a Set of Trees
#'
#' @description
#' This is a wrapper for TransPhylo::infer_multittree_share_param. 
#' Can use it to run TransPhylo on a set of trees, either all MCC trees
#' with parameter sharing, or a sample of posterior BEAST2 trees without
#' parameter sharing (i.e., a sample of BEAST2 posterior trees for a single
#' cluster, to incorporate some phylogenetic uncertainty).
#'
#' @param trees An `ape` multiPhylo object containing the trees to be
#' analyzed. If `type` is "mcctrees", this should be all MCC trees for
#' all clusters, and the names of the trees need to match the names
#' of the trees in the cluster dictionary. If `type` is "trees_sample",
#' this should be a sample of posterior BEAST2 trees for a single cluster,
#' and the cluster name should be provided as an argument.
#' @param cluster_name If type is "trees_sample", this should be the name
#' of the cluster for which the trees are being analyzed. This should match
#' the cluster name in the cluster dictionary. If type is "mcctrees", this
#' argument is ignored.
#' @param type A character string specifying the type of trees to be
#' analyzed. Can be "mcctrees" or "trees_sample". If "mcctrees", this
#' function will run TransPhylo on all MCC trees for all clusters. If
#' "trees_sample", this function will run TransPhylo on a sample of
#' posterior BEAST2 trees for a single cluster.
#' @param cluster_dict A data frame containing `sample_id`, `cluster_name`,
#' and `collectdt` for each sequence, as output by `assign_snp_clusters`.
#' @param out_dir A character string specifying the directory where
#' the TransPhylo results will be saved.
#' @param output_name A character string specifying the name of the output
#' file.
#' @param ... Additional arguments to be passed to
#' `TransPhylo::infer_multittree_share_param`.
#'
#' @returns A data frame containing the probability of each sample being an
#' infection source.
#'
#' @import ape
#' @import TransPhylo
#' @import lubridate
#'
#' @export

run_TransPhylo <- function(trees,
                           cluster_name = NULL,
                           type = c("mcctrees", "trees_sample"),
                           cluster_dict,
                           out_dir = "TransPhylo",
                           output_name = "tp_res",
                           ...) {

  # Prep data for TP run
  args <- prep_for_TransPhylo(trees = trees,
                              cluster_name = cluster_name,
                              type = type,
                              cluster_dict = cluster_dict,
                              out_dir = out_dir,
                              output_name = output_name,
                              ...)

  # Run TransPhylo multitree
  tp_res <- do.call(TransPhylo::infer_multittree_share_param, args)

  # Save resTransPhylo object(s)
  if (!dir.exists(out_dir)) dir.create(out_dir)
  # For multitree, TransPhylo forgets that it's supposed to be a
  # list of resTransPhylo objects
  tp_res <- lapply(tp_res, function(res) {
    class(res) <- "resTransPhylo"
    res
  })
  saveRDS(tp_res, file = file.path(out_dir,
                                   paste0(output_name, ".rds")))

  # Get prob_source for all samples
  prob_source <- get_prob_source(tp_res, type = type)

  # And save prob_source
  output_name <- paste0(output_name, "_prob_source.rds")
  saveRDS(prob_source, file = file.path(out_dir, output_name))
  return(prob_source)
}

#' @title Prepare Data for TransPhylo
#'
#' @inheritParams run_TransPhylo
#'
#' @returns Arguments to be passed into
#' TransPhylo::infer_multittree_share_param().

prep_for_TransPhylo <- function(trees,
                                cluster_name = NULL,
                                type = c("mcctrees", "trees_sample"),
                                cluster_dict,
                                out_dir = "TransPhylo",
                                output_name = "tp_res",
                                ...) {
    # Some TP arguments are set in this function - ignore
  # if they are passed in
  if ("ptree_lst" %in% names(list(...))) {
    warning("Ignoring user-provided input ptree_list, which is set in this function.")
  }
  if ("dateT" %in% names(list(...))) {
    warning("Ignoring user-provided input dateT, which is set in this function.")
  }

  # Ensure names of trees are in provided cluster dictionary
  # (Not all clusters in the dictionary need a tree, but all
  # trees need to be in the dictionary)
  if (type == "mcctrees") {
    if (!all(names(trees) %in% cluster_dict$cluster_name)) {
      stop("Not all tree names are in the cluster dictionary.")
    }
  } else if (type == "trees_sample") {
    if (is.null(cluster_name)) {
      stop("Cluster name must be provided for type 'trees_sample'.")
    }
    if (!cluster_name %in% cluster_dict$cluster_name) {
      stop("Cluster name not in cluster dictionary.")
    }
  }

  # Convert dates in the cluster dict to decimal dates
  cluster_dict$collectdt <- as.Date(cluster_dict$collectdt)
  cluster_dict$collectdt <- lubridate::decimal_date(cluster_dict$collectdt)

  # Convert phylo (ape) -> ptree (TransPhylo)
  if (type == "mcctrees") {
    # Get the date of the last sample for each cluster
    dates_last_sample <- aggregate(collectdt ~ cluster_name, data = cluster_dict, max)
    # Reorder dates last sample to match tree order
    dates_last_sample <- dates_last_sample[match(names(trees),
                                                 dates_last_sample$cluster_name),
                                           "collectdt"]
    # Convert ape phylo trees to TP ptrees using date_last_sample
    ptree_list <- mapply(FUN = function(tree, date_last_sample) {
      TransPhylo::ptreeFromPhylo(tree, date_last_sample)
    }, trees, dates_last_sample, SIMPLIFY = FALSE)
    # Get the latest date for the TransPhylo run
    date_last_sample <- max(dates_last_sample)
  } else if (type == "trees_sample") {
    # Pull the correct cluster from the cluster dictionary
    # & get date of last sample
    date_last_sample <- cluster_dict[cluster_dict$cluster_name == cluster_name, ]
    date_last_sample <- max(date_last_sample$collectdt)
    # Convert to ptrees
    ptree_list <- lapply(trees, FUN = function(tree) {
      TransPhylo::ptreeFromPhylo(tree, date_last_sample)
    })
  }
  # Force no sharing for trees_sample
  args <- list(...)
  if (type == "trees_sample") {
    args$share <- NULL
  }
  args$ptree_lst <- ptree_list
  args$dateT <- date_last_sample + 0.01
  return(args)
}

#' @title Get Infector Probabilities from TransPhylo Results
#'
#' @param tp_res A list of resTransPhylo objects.
#' @param type A character string specifying the type of trees that were
#' analyzed. Can be "mcctrees" or "trees_sample".
#'
#' @returns A data frame containing the probability of each sample being an
#' infection source.

get_prob_source <- function(tp_res,
                            type = c("mcctrees", "trees_sample")) {
    # Get prob source for all samples & save
  names <- lapply(tp_res, function(res) {
    res[[1]]$ctree$nam
  })
  offspring <- mapply(FUN = function(res, names) {
    TransPhylo::getOffspringDist(res, names, burnin = .5)
  }, tp_res, names, SIMPLIFY = FALSE)
  prob_source <- lapply(offspring, function(offspring) {
    rowSums(offspring != 0) / dim(offspring)[2]
  })
  prob_source <- mapply(FUN = function(prob_source, names) {
    data.frame(names = names, prob_source = prob_source,
               stringsAsFactors = FALSE)
  }, prob_source, names, SIMPLIFY = FALSE)
  prob_source <- do.call(rbind, prob_source)
  if (type == "trees_sample") {
    prob_source <- aggregate(prob_source ~ names,
                             data = prob_source, FUN = mean)
  }
  colnames(prob_source) <- c("SampleID", "prob_source")
  return(prob_source)
}