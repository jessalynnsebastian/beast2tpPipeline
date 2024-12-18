# beast2tpPipeline Package

## Introduction

This package was created to work with data from the Kopanyo study, namely fasta files of TB sequences and
associated metadata. It defines a clear pipeline from fasta files of sequences and associated metadata, through
creating SNP clusters, running BEAST2, running TransPhylo, and doing regression on the results.

In [Using genetic data to identify transmission risk factors: Statistical assessment and application to tuberculosis transmission](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010696),
Goldstein et al. show that this pipeline is very biased; the BEAST2-TransPhylo pipeline has low sensitivity of identifying infection sources and biases regression coefficients toward zero. This package was built to facilitate
exploration of the qualities of this pipeline.

This package is built such that each step in the pipeline has an associated function. The steps and functions are:

1. **Assigning SNP clusters**: `assign_snp_clusters` takes a set of sequences and metadata, and assigns sequences to SNP clusters.
2. **Creating BEAST2 clusters**: `create_BEAST2_clusters` takes the SNP clusters and separates the sequence data by cluster, keeping only the SNPs if desired.
3. **Creating XML files for BEAST2**: `create_cluster_xml` takes the SNP data and creates an XML file for BEAST2 to run.
4. **Running BEAST2**: `run_beast2` runs BEAST2 on the XML file. This is where we take the DNA sequences and infer a timed phylogenetic tree, which tells us about the ancestry of the pathogen.
5. **Getting the MCC trees**: `get_mcctree` takes the posterior trees from BEAST2 and creates a maximum clade credibility tree - a summary tree, like a mean or median.
(Alternatively, **Getting a Sample of Posterior Trees**: `sample_BEAST2_trees` takes the posterior trees from BEAST2
and subsamples a specified number of them.)
6. **Running TransPhylo on the MCC trees**: `run_TransPhylo` runs TransPhylo on the MCC trees, sharing some parameters that can be simultaneously estimated. This is where we take the timed phylogeny and use it to try to infer who-infected-whom.
(Alternatively, you can use `run_TransPhylo` to run TransPhylo on a subsample of BEAST2 trees, without parameter sharing.)
7. **Running a regression on the TransPhylo results**: `regression` runs linear or logistic regression using the TransPhylo results. This is where we draw inference about the relationship between some covariates and individual transmission probabilities.

## Getting Started

Install the `beast2tpPipeline` package from GitHub:

```r
devtools::install_github("jessalynnsebastian/beast2tpPipeline", build_vignettes = TRUE)
```

There are some dependencies that should be automatically downloaded. Then, load the package:

```r
library(beast2tpPipeline)
```

You will also need to have BEAST2 installed on your computer. You can download it [here](https://www.beast2.org/). If you use
a package manager like homebrew, you can install BEAST2 with `brew install beast2` in your terminal. This will install a suite 
of programs, including BEAST2, BEAUti, and TreeAnnotator. You will probably also want to install Tracer to look at traceplots
and other diagnostics of your BEAST2 runs. You can download Tracer [here](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).
