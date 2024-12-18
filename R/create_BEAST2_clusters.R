#' @title Create Clusters for BEAST2 Input
#'
#' @description
#' Using SNP cluster information and DNA sequences, separate sequences into
#' separate DNAbin objcts by cluster. If desired, keep only SNPs.
#'
#' @param seqs A matrix of DNA sequences, where rownames correspond
#' to sample IDs, or an `ape::DNAbin` object.
#' @param cluster_assignments A data frame of cluster assignments, including
#' at least `sample_id` and `cluster_name`, as output by `assign_snp_clusters`.
#' Must also include `collectdt` if writing a cluster dictionary.
#' @param min_cluster_size The number of samples required to keep a cluster.
#' If the cluster is too small, discard it.
#' @param min_varsites The number of variable sites required to keep a cluster.
#' If there aren't enough variable sites, discard the cluster.
#' @param snps_only If TRUE, keep only SNPs in the output DNAbin objects.
#' @param constant_sites A string of constant sites to add to the beginning of
#' each sequence. This is necessary for BEAST2 input if using only SNPs. Defaults
#' to an empty string.
#' @param cluster_dictionary_file Optionally, a file path to write a cluster
#' dictionary to, containing cluster names, sizes, and varsites.
#' @param fasta_dir Optionally, a directory to write FASTA files to for each
#' cluster. FASTA files are automatically named after the cluster names.
#'
#' @returns A named list of matrices containing SNPs.
#'
#' @import ape
#' @import stats
#' @import utils
#'
#' @export

create_BEAST2_clusters <- function(seqs,
                                   cluster_assignments,
                                   min_cluster_size = 4,
                                   min_varsites = 8,
                                   snps_only = TRUE,
                                   constant_sites = "",
                                   cluster_dictionary_file = NULL,
                                   fasta_dir = NULL) {
  # Check that sequences are in correct format
  if (!is.matrix(seqs)) {
    if (inherits(seqs, "DNAbin")) {
      seqs <- ape::as.character.DNAbin(seqs)
    } else {
      stop("Sequences not in correct format. Make sure you are reading
          sequences in with ape::read.dna([FILE], format = 'fasta',
          as.character = TRUE).")
    }
  }

  # Cut clusters that are too small
  cluster_sizes <- table(cluster_assignments$cluster_name)
  small_clusters <- names(cluster_sizes)[cluster_sizes < min_cluster_size]
  cluster_assignments$cluster_name[cluster_assignments$cluster_name %in% small_clusters] <- NA
  cluster_assignments <- na.omit(cluster_assignments)

  # Match cluster assignments to sequences via sample ID
  seq_cluster_df <- merge(cluster_assignments, seqs,
                          by.x = "sample_id", by.y = "row.names")
  cluster_dna_list <- split(
    seq_cluster_df,
    seq_cluster_df$cluster_name
  )
  names(cluster_dna_list) <- lapply(cluster_dna_list, function(cluster) {
    unique(cluster$cluster_name)
  })

  # Now that clusters are split, remove unnecessary columns & leave only seqs
  cluster_dna_list <- lapply(cluster_dna_list, function(cluster) {
    rownames(cluster) <- cluster$sample_id
    cluster <- cluster[, -c(1, 2, 3), drop = FALSE]
    cluster
  })

  if (snps_only) {
    # Keep only SNPs using seg.sites
    cluster_dna_list <- lapply(cluster_dna_list, function(cluster_df) {
      cluster_df <- ape::as.DNAbin(as.matrix(cluster_df))
      cluster_df <- cluster_df[, ape::seg.sites(cluster_df)]
      if (dim(cluster_df)[2] < min_varsites) {
        return(NULL)
      }
      # Add const sites for BEAST2
      cluster_df <- ape::as.alignment(cluster_df)
      cluster_df$seq <- sapply(cluster_df$seq, function(x) {
        paste0(constant_sites, x)
      })
      cluster_df <- ape::as.DNAbin(cluster_df)
      ape::as.character.DNAbin(cluster_df)
    })
  }

  # Remove NULL clusters
  cluster_dna_list <- cluster_dna_list[!sapply(cluster_dna_list, is.null)]

  # Create a dictionary of cluster information
  if (!is.null(cluster_dictionary_file)) {
    if (is.null(cluster_assignments$collectdt)) {
      stop("Please provide collection dates to write a cluster dictionary.")
    }
    cluster_info <- data.frame(cluster_name = names(cluster_dna_list),
                               cluster_size = lapply(cluster_dna_list, length),
                               cluster_snps = lapply(cluster_dna_list, function(cluster) {
                                 length(cluster[[1]])
                               }),
                               cluster_sampleids = lapply(cluster_dna_list, function(cluster) {
                                 rownames(cluster)
                               }),
                               date_last_sample = lapply(names(cluster_dna_list), function(cluster_name) {
                                 collectdts <- cluster_assignments$collectdt[cluster_assignments$cluster_name == cluster_name]
                                 max(as.Date(collectdts))
                               }))
    write.csv(cluster_info, cluster_dictionary_file, row.names = FALSE)
  }

  # Write FASTA files for each cluster
  if (!is.null(fasta_dir)) {
    for (cluster_name in names(cluster_dna_list)) {
      fasta_path <- file.path(fasta_dir, paste0(cluster_name, ".fasta"))
      ape::write.dna(cluster_dna_list[[cluster_name]], file = fasta_path, format = "fasta")
    }
  }
  
  return(cluster_dna_list)
}