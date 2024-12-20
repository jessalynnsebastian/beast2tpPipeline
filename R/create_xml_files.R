#' @title Create BEAST2 XML Files from a Template
#'
#' @description
#' Use a template XML file to create new XML files for BEAST2 input.
#' This function can replace sequences, sampling dates, mcmc iterations,
#' frequency of storing trees, and, for a uniform clock rate, the minimum
#' and maximum clock rates. 
#' 
#' This was designed just to be used with the provided template file, though
#' it could theoretically be used with others. The template file was written
#' for tuberculosis data. It contains an adjustment for ascertainment bias due
#' to the SNPs, a Gamma site model, and HKY substitution model, strict clock,
#' Coalescent constant population model, uniform clock rate prior (min/max editable),
#' and lognormal(1, 1.25) freqParameter, kappa, and popsize priors.
#' 
#' The template xml contains placeholders:
#' - SNP_FILE_NAME_HERE for the name of the SNP file
#' - ALIGNMENT_INFORMATION_HERE for the sequence information
#' - DATE_INFORMATION_HERE for the sampling dates
#' - CLOCKRATE_MINIMUM_HERE for the minimum uniform clock rate
#' - CLOCKRATE_MAXIMUM_HERE for the maximum uniform clock rate 
#' - CLOCKRATE_INITIAL_HERE for the initial clock rate
#' - MCMC_ITERATIONS_HERE for the number of MCMC iterations
#' - STORE_EVERY_HERE for the frequency of storing trees
#' 
#' 
#' @param path_to_template The path to the template XML file.
#' @param seqs_list A list of matrices of DNA sequences (or SNPs), where
#' rownames correspond to sample IDs.
#' @param cluster_assignments A data frame of cluster assignments, including
#' `sample_id`, `cluster_name`, and `collectdt`, as output by
#' `assign_snp_clusters`. All sample IDs in the `seqs` matrices must be
#' present in this data frame.
#' @param mcmc_iterations The number of MCMC iterations to run.
#' @param store_every The frequency of storing trees.
#' @param min_clockrate The minimum clock rate for the uniform clock
#' rate prior. Defaults to tuberculosis genome clock rate. DO NOT adjust
#' for SNPs here.
#' @param max_clockrate The maximum clock rate for the uniform clock
#' rate prior. Defaults to tuberculosis genome clock rate.
#' @param init_clockrate The initial clock rate for the clock rate.
#' @param whole_genome_length The length of the whole genome. Used
#' to adjust the clock rate for SNPs. Defaults to TB length.
#' @param out_dir The directory to write the output XML file. Files
#' are automatically named after their FASTA names.
#'
#' @returns NULL
#'
#' @import ape
#'
#' @export
create_xml_files <- function(path_to_template = system.file("xml_template/0_xml_template.xml",
                                                                package = "beast2tpPipeline"),
                             seqs_list,
                             cluster_assignments,
                             mcmc_iterations = 20000000,
                             store_every = 5000,
                             min_clockrate = 10^(-8),
                             max_clockrate = 5 * 10^(-7),
                             init_clockrate = 10 * min_clockrate,
                             whole_genome_length = 4.2 * 10^6,
                             out_dir) {
  # Get xmls for all clusters
  if (length(seqs_list) == 0) {
    message("No sequences provided.")
    return(NULL)
  }
  cluster_assignments_list <- split(cluster_assignments,
                                    cluster_assignments$cluster_name)
  cluster_names <- names(seqs_list)

  # Get init_clockrate
  if (missing(init_clockrate)) {
    init_clockrate <- 10 * min_clockrate
  }

  for (name in cluster_names) {
    message(paste("Creating XML for", name))
    create_cluster_xml(path_to_template,
                       seqs_list[[name]],
                       name,
                       cluster_assignments_list[[name]],
                       mcmc_iterations,
                       store_every,
                       min_clockrate,
                       max_clockrate,
                       init_clockrate,
                       whole_genome_length,
                       out_dir)
    message(paste("XML created for", name, "in", out_dir))
  }
}

#' @title Create a BEAST2 XML File from a Template
#' 
#' @param cluster_seqs A matrix of DNA sequences (or SNPs), where
#' rownames correspond to sample IDs.
#' @param cluster_name A character string specifying the name of the cluster.
#' @inheritParams create_xml_files
#'
#' @returns NULL
#'
#' @import ape

create_cluster_xml <- function(path_to_template,
                               cluster_seqs,
                               cluster_name,
                               cluster_assignments,
                               mcmc_iterations,
                               store_every,
                               min_clockrate,
                               max_clockrate,
                               init_clockrate,
                               whole_genome_length,
                               out_dir) {
  if (is.null(out_dir)) {
    stop("Please provide an output directory.")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  if (!file.exists(path_to_template)) {
    stop("Template XML file not found.")
  }
  if (is.null(cluster_seqs)) {
    stop("Please provide SNP data.")
  }
  if (is.null(cluster_name)) {
    stop("Please provide a cluster name.")
  }
  if (is.null(cluster_assignments)) {
    stop("Please provide cluster assignments. You can obtain them using `assign_snp_clusters`.")
  }
  # Read in the template XML file & replace name of SNP file
  xml_template <- readLines(path_to_template)
  new_xml <- gsub("SNP_FILE_NAME_HERE", cluster_name, xml_template)

  # Replace the SNP data
  if (!is.matrix(cluster_seqs)) {
    if (inherits(cluster_seqs, "DNAbin")) {
      seqs_matrix <- ape::as.character.DNAbin(cluster_seqs)
    } else {
      stop("seqs not in correct format. Make sure you are reading
          sequences in with ape::read.dna([FILE], format = 'fasta',
          as.character = TRUE).")
    }
  } else {
    seqs_matrix <- cluster_seqs
  }
  snp_length <- ncol(seqs_matrix) - 4 # subtract 4 for the constant sites
  sequence_string <- '<sequence id="seq_SEQUENCE_NAME" spec="Sequence" taxon="SEQUENCE_NAME" totalcount="4" value="SNP_SEQUENCE"/>'
  # Note that SNP_SEQUENCE needs to contain the 4 const sites, acgt
  for (i in seq_len(nrow(seqs_matrix))) {
    insert_string <- gsub("SEQUENCE_NAME",
                          rownames(seqs_matrix)[i],
                          sequence_string)
    insert_string <- gsub("SNP_SEQUENCE",
                          paste(seqs_matrix[i, ], collapse = ""),
                          insert_string)
    if (i != nrow(seqs_matrix)) {
      insert_string <- paste(insert_string,
                            "\n\t\tALIGNMENT_INFORMATION_HERE",
                            sep = "")
    }
    new_xml <- gsub("ALIGNMENT_INFORMATION_HERE",
                    insert_string,
                    new_xml,
                    fixed = TRUE)
  }

  # Replace the sampling dates
  cluster_assignments <- cluster_assignments[cluster_assignments$cluster_name == cluster_name, ]
  sampling_dates <- setNames(cluster_assignments$collectdt,
                            cluster_assignments$sample_id)
  sampling_dates <- paste(names(sampling_dates), sampling_dates, sep = "=")
  new_xml <- gsub("DATE_INFORMATION_HERE",
                  paste(sampling_dates, collapse = ","),
                  new_xml)

  # Replace the clock rate info
  min_unif_clockrate <- (min_clockrate * whole_genome_length) / snp_length
  max_unif_clockrate <- (max_clockrate * whole_genome_length) / snp_length
  init_clockrate <- (init_clockrate * whole_genome_length) / snp_length

  new_xml <- gsub("CLOCKRATE_MINIMUM_HERE", min_unif_clockrate, new_xml)
  new_xml <- gsub("CLOCKRATE_MAXIMUM_HERE", max_unif_clockrate, new_xml)
  new_xml <- gsub("CLOCKRATE_INITIAL_HERE", init_clockrate, new_xml)

  # Change mcmc iterations and store every
  new_xml <- gsub("MCMC_ITERATIONS_HERE", mcmc_iterations, new_xml)
  new_xml <- gsub("STORE_EVERY_HERE", store_every, new_xml)

  # Write the new XML file
  out_path <- file.path(out_dir, paste0(cluster_name, ".xml"))
  writeLines(new_xml, out_path)
}