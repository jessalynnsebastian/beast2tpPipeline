#' @title Assign SNP Clusters to Samples
#'
#' @description
#' This function assigns SNP clusters to samples based on the SNP distance and
#' some threshold, or, alternatively, based on transmission clusters as in the
#' Stimson et al. paper, \href{https://academic.oup.com/mbe/article/36/3/587/5300248}{Beyond the SNP Threshold: Identifying Outbreak Clusters Using Inferred Transmissions}.
#' 
#' @param seqs Matrices containing sequences.
#' @param collectdts A numeric vector of collection dates, with names corresponding
#' to the rownames of the `seqs` matrix.
#' @param snp_matrix If SNP distances have been previously computed, include
#'  this argument to avoid re-computing.
#' @param threshold A numeric value of the SNP threshold, or the transmission
#'  threshold if using transmission clusters (see `transcluster` for more info).
#' @param clockrate A numeric value of the clock rate if using transmission
#'  clusters.
#' @param transm_rate A numeric value of the transmission rate if using
#'  transmission clusters.
#'
#' @returns A data frame containing sample ID, cluster name, and collection date
#' for the sample.
#'
#' @import transcluster
#' @import ape
#'
#' @export

assign_snp_clusters <- function(seqs = NULL,
                                collectdts,
                                snp_matrix = NULL,
                                threshold = 5,
                                clockrate = -1,
                                transm_rate = -1) {
  # Get SNP distances
  if (is.null(snp_matrix)) {
    snp_matrix <- as.matrix(ape::dist.gene(seqs))
  }
  # Split into clusters
  seq_sampleids <- rownames(snp_matrix)
  collectdts <- collectdts[seq_sampleids]
  collectdts_yyyy_mm_dd <- as.Date(collectdts)
  collectdts <- lubridate::decimal_date(collectdts_yyyy_mm_dd)
  clust_model <- transcluster::createModel(ids = seq_sampleids,
                                           dates = collectdts,
                                           snpMatrix = snp_matrix)
  if (clockrate == -1 || transm_rate == -1) {
    ## SNP-based clusters
    message("Creating SNP-based clusters")
    clust_model <- transcluster::setSNPThresholds(clust_model, threshold)
    clusters <- transcluster::makeSNPClusters(clust_model, writeFile = FALSE)
  } else {
    ## Transmission-based clusters
    message("Creating transmission-based clusters")
    clust_model <- transcluster::setParams(clust_model, lambda = clockrate,
                                           beta = transm_rate)
    clust_model <- transcluster::setTransThresholds(clust_model,
                                                    threshold)
    clusters <- transcluster::makeTransClusters(clust_model, writeFile = FALSE)
  }
  # Label sequences with their cluster
  clust_labels <- data.frame(sample_id = seq_sampleids,
                             cluster_name = paste0("cluster", clusters[[1]]),
                             collectdt = collectdts_yyyy_mm_dd,
                             row.names = NULL)
  return(clust_labels)
}
