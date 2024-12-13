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
#' @param gentime_shape A numeric specifying the shape parameter of the gamma
#' distribution for the generation time.
#' @param gentime_scale A numeric specifying the scale parameter of the gamma
#' distribution for the generation time.
#' @param sampling_shape A numeric specifying the shape parameter of the gamma
#' distribution for the sampling time. Default is the same as `gentime_shape`.
#' @param sampling_scale A numeric specifying the scale parameter of the gamma
#' distribution for the sampling time. Default is the same as `gentime_scale`.
#' @param prior_sampfrac_a A numeric specifying the first shape parameter of the beta
#' distribution for the sampling fraction.
#' @param prior_sampfrac_b A numeric specifying the second shape parameter of the beta
#' distribution for the sampling fraction.
#' @param start_off_p A numeric specifying the starting value for the off.p parameter.
#' @param start_neg A numeric specifying the starting value for the neg parameter.
#' @param start_pi A numeric specifying the starting value for the pi parameter.
#' @param start_off_r A numeric specifying the starting value for the off.r parameter.
#' @param share A character vector specifying which parameters to share across
#' clusters. Defaults to "neg" and "off.r". If `type` is "trees_sample", this
#' parameter is ignored.
#' @param update_neg A logical specifying whether to update the neg parameter.
#' @param update_off_r A logical specifying whether to update the off.r parameter.
#' @param update_off_p A logical specifying whether to update the off.p parameter.
#' @param update_pi A logical specifying whether to update the pi parameter.
#' @param optiStart A numeric specifying type of optimization to apply to
#' MCMC start point.
#' @param verbose A logical specifying whether to print verbose output.
#' @param mcmc_iterations A numeric specifying the number of MCMC iterations.
#' @param thinning A numeric specifying the thinning parameter for the MCMC.
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
                           gentime_shape = 10,
                           gentime_scale = 1 / 10,
                           sampling_shape = 10,
                           sampling_scale = 1/10,
                           prior_sampfrac_a = 1,
                           prior_sampfrac_b = 19,
                           start_off_p = .5,
                           start_neg = 1.48,
                           start_pi = .5,
                           start_off_r = 1,
                           share = c("neg", "off.r"),
                           update_neg = TRUE,
                           update_off_r = TRUE,
                           update_off_p = FALSE,
                           update_pi = TRUE,
                           optiStart = 2,
                           verbose = FALSE,
                           mcmc_iterations = 100000,
                           thinning = 10) {
  # Force no sharing for trees_sample
  if (type == "trees_sample") share <- NULL

  # Ensure names of trees are in provided cluster dictionary
  # (Not all clusters in the dictionary need a tree, but all
  # trees need to be in the dictionary)
  if (type == "mcctrees") {
    if (!all(names(trees) %in% cluster_dict$cluster_name)) {
      stop("Not all tree names are in the cluster dictionary.")
    }
  } else if (type == "trees_sample") {
    if (!cluster_name %in% cluster_dict$cluster_name) {
      stop("Cluster name not in cluster dictionary.")
    }
  }

  # Convert phylo (ape) -> ptree (TransPhylo)
  if (type == "mcctrees") {
    # Get the date of the last sample for each cluster
    dates_last_sample <- aggregate(collectdt ~ cluster_name, data = cluster_assignments, max)
    # Reorder dates last sample to match tree order
    dates_last_sample <- dates_last_sample[match(names(trees),
                                                 dates_last_sample$cluster_name),
                                          "collectdt"]
    # Convert ape phylo trees to TP ptrees using date_last_sample
    ptree_list <- mapply(FUN = function(tree, date_last_sample) {
      TransPhylo::ptreeFromPhylo(tree, date_last_sample)
    }, trees, dates_last_sample, SIMPLIFY = FALSE)
  } else if (type == "trees_sample") {
    # Pull the correct cluster from the cluster dictionary
    # & get date of last sample
    cluster_dict <- cluster_dict[cluster_dict$cluster_name == cluster_name, ]
    # Convert to ptrees
    ptree_list <- lapply(trees, FUN = function(tree) {
      TransPhylo::ptreeFromPhylo(tree, date_last_sample)
    })
  }
  # Get the latest date from the cluster dictionary for the
  # TransPhylo run
  date_last_sample <- max(cluster_dict$collectdt)
  # Run TransPhylo multitree
  tp_res <- TransPhylo::infer_multittree_share_param(ptree_list,
                                                     mcmcIterations = mcmc_iterations,
                                                     thinning = thinning,
                                                     w.shape = gentime_shape,
                                                     w.scale = gentime_scale,
                                                     ws.shape = sampling_shape,
                                                     ws.scale = sampling_scale,
                                                     startOff.r = start_off_r,
                                                     startOff.p = start_off_p,
                                                     startNeg = start_neg,
                                                     startPi = start_pi,
                                                     prior_pi_a = prior_sampfrac_a,
                                                     prior_pi_b = prior_sampfrac_b,
                                                     updateNeg = update_neg,
                                                     updateOff.r = update_off_r,
                                                     updateOff.p = update_off_p,
                                                     updatePi = update_pi,
                                                     share = share,
                                                     optiStart = optiStart,
                                                     dateT = date_last_sample + .01,
                                                     verbose = verbose)

  # Save resTransPhylo object(s)
  if(!dir.exists(output_directory)) dir.create(output_directory)
  saveRDS(tp_res, file = file.path(output_directory,
                                   paste0(output_name, ".rds")))

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
  # And save the vector
  output_name <- paste0(output_name, "_prob_source.rds")
  saveRDS(prob_source, file = file.path(output_directory, output_name))
}