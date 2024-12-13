#' @title Create BEAST2 XML File
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
#' @param seqs A matrix of DNA sequences (or SNPs), where rownames correspond
#' to sample IDs, or an `ape::DNAbin` object.
#' @param cluster_name The name of the cluster.
#' @param sampling_dates A named vector of sampling dates, where
#' the names are the corresponding sample IDs. Dates should be in
#' either decimal years or in the format "YYYY-MM-DD".
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

create_cluster_xml <- function(path_to_template = system.file("xml_template/0_xml_template.xml",
                                                              package = "beast2tpPipeline"),
                               snps,
                               cluster_name,
                               sampling_dates,
                               mcmc_iterations = 10000000,
                               store_every = 5000,
                               min_clockrate = 10^(-8),
                               max_clockrate = 5 * 10^(-7),
                               init_clockrate = 10 * min_clockrate,
                               whole_genome_length = 4.2 * 10^6,
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
  if (is.null(snps)) {
    stop("Please provide SNP data.")
  }
  if (is.null(cluster_name)) {
    stop("Please provide a cluster name.")
  }
  if (is.null(sampling_dates)) {
    stop("Please provide sampling dates. Make sure that they are a named vector, where names are sample IDs.")
  }
  # Read in the template XML file & replace name of SNP file
  xml_template <- readLines(path_to_template)
  new_xml <- gsub("SNP_FILE_NAME_HERE", cluster_name, xml_template)

  # Replace the SNP data
  if (!is.matrix(snps)) {
    if (class(snps) == "DNAbin") {
      snps_matrix <- ape::as.character.DNAbin(snps)
    } else {
      stop("SNPs not in correct format. Make sure you are reading
          sequences in with ape::read.dna([FILE], format = 'fasta',
          as.character = TRUE).")
    }
  } else {
    snps_matrix <- snps
  }
  snp_length <- ncol(snps_matrix) - 4 # subtract 4 for the constant sites
  sequence_string <- '<sequence id="seq_SEQUENCE_NAME" spec="Sequence" taxon="SEQUENCE_NAME" totalcount="4" value="SNP_SEQUENCE"/>'
  # Note that SNP_SEQUENCE needs to contain the 4 const sites, acgt
  for (i in seq_len(nrow(snps_matrix))) {
    insert_string <- gsub("SEQUENCE_NAME",
                          rownames(snps_matrix)[i],
                          sequence_string)
    insert_string <- gsub("SNP_SEQUENCE",
                          paste(snps_matrix[i, ], collapse = ""),
                          insert_string)
    if (i != nrow(snps_matrix)) {
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
  sampling_dates <- paste(names(sampling_dates), sampling_dates, sep = "=")
  new_xml <- gsub("DATE_INFORMATION_HERE",
                  paste(sampling_dates, collapse = ","),
                  new_xml)

  # Replace the clock rate info
  min_unif_clockrate <- (min_clockrate * whole_genome_length) / snp_length
  max_unif_clockrate <- (max_clockrate * whole_genome_length) / snp_length

  new_xml <- gsub("CLOCKRATE_MINIMUM_HERE", min_unif_clockrate, new_xml)
  new_xml <- gsub("CLOCKRATE_MAXIMUM_HERE", max_unif_clockrate, new_xml)
  if (missing(init_clockrate)) {
    init_clockrate <- 10 * min_unif_clockrate
  }
  new_xml <- gsub("CLOCKRATE_INITIAL_HERE", init_clockrate, new_xml)

  # Change mcmc iterations and store every
  new_xml <- gsub("MCMC_ITERATIONS_HERE", mcmc_iterations, new_xml)
  new_xml <- gsub("STORE_EVERY_HERE", store_every, new_xml)

  # Write the new XML file
  out_path <- file.path(out_dir, paste0(cluster_name, ".xml"))
  writeLines(new_xml, out_path)
}
