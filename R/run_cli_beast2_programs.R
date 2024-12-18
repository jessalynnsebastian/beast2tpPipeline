# Functions that call command-line BEAST2 programs.
# For more info on running BEAST2 via command line, see:
# https://www.beast2.org/2019/09/26/command-line-tricks.html

#' @title Run Command-Line BEAST2
#'
#' @description
#' Run BEAST2 from an XML file. Tree files, log files, etc. will be written
#' to the same directory as the XML file.
#'
#' @param input_xml_path Character string. Path to the BEAST2 XML file.
#' @param beast2_path Character string. Path to the BEAST2 executable.
#'
#' @returns NULL
#'
#' @export
run_beast2 <- function(input_xml_path,
                       beast2_path = '/Applications/"BEAST 2.7.7"/bin/beast') {
  opts <- "-working -overwrite"
  # changes output dir to working directory
  # and sets to overwrite existing log files
  cmd <- paste(beast2_path, opts, input_xml_path)
  system(cmd)
}

#' @title Get Maximum Clade Credibility Trees Using Command-Line TreeAnnotator
#'
#' @description
#' Run TreeAnnotator to get the maximum clade credibility tree from a BEAST2
#' tree file.
#'
#' @param input_treesfile Character string. Path to the BEAST2 .trees file.
#' @param output_dir Character string. Path to the desired output directory
#' for the MCC tree nexus file.
#' @param beast_iterations Numeric. Number of iterations in the BEAST2 run.
#' @param burnin_fraction Numeric. Fraction of the MCMC chain to discard as
#' burn-in.
#' @param heights Character string. Node heights for the tree. Beware that if
#' using anything other than common ancestor, could get some weird trees that
#' won't work with TransPhylo.
#' @param treeannotator_path Character string. Path to the TreeAnnotator
#' executable.
#'
#' @returns NULL
#'
#' @export
get_mcctree <- function(input_treesfile, output_dir,
                        beast_iterations = 20000000,
                        burnin_fraction = 1 / 2,
                        heights = "CA",
                        treeannotator_path = '/Applications/"BEAST 2.7.7"/bin/treeannotator') {
  # Get command as a character string
  out_name <- gsub(".trees", ".nexus", basename(input_treesfile))
  out_name <- gsub("^[^-]+", "mcctree", out_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  cmd <- paste(treeannotator_path, "-burnin", 100 * burnin_fraction,
               "-height", heights, input_treesfile, file.path(output_dir,
                                                              out_name))
  # Run command-line treeannotator
  system(cmd)
}