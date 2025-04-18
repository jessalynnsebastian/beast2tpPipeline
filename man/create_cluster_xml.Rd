% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_xml_files.R
\name{create_cluster_xml}
\alias{create_cluster_xml}
\title{Create a BEAST2 XML File from a Template}
\usage{
create_cluster_xml(
  path_to_template,
  cluster_seqs,
  cluster_name,
  cluster_assignments,
  mcmc_iterations,
  store_every,
  min_clockrate,
  max_clockrate,
  init_clockrate,
  whole_genome_length,
  out_dir
)
}
\arguments{
\item{path_to_template}{The path to the template XML file.}

\item{cluster_seqs}{A matrix of DNA sequences (or SNPs), where
rownames correspond to sample IDs.}

\item{cluster_name}{A character string specifying the name of the cluster.}

\item{cluster_assignments}{A data frame of cluster assignments, including
\code{sample_id}, \code{cluster_name}, and \code{collectdt}, as output by
\code{assign_snp_clusters}. All sample IDs in the \code{seqs} matrices must be
present in this data frame.}

\item{mcmc_iterations}{The number of MCMC iterations to run.}

\item{store_every}{The frequency of storing trees.}

\item{min_clockrate}{The minimum clock rate for the uniform clock
rate prior. Defaults to tuberculosis genome clock rate. DO NOT adjust
for SNPs here.}

\item{max_clockrate}{The maximum clock rate for the uniform clock
rate prior. Defaults to tuberculosis genome clock rate.}

\item{init_clockrate}{The initial clock rate for the clock rate.}

\item{whole_genome_length}{The length of the whole genome. Used
to adjust the clock rate for SNPs. Defaults to TB length.}

\item{out_dir}{The directory to write the output XML file. Files
are automatically named after their FASTA names.}
}
\description{
Create a BEAST2 XML File from a Template
}
\keyword{internal}
