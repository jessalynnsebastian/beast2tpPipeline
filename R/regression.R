#' @title Run Linear or Logistic Regression with TransPhylo Results
#'
#' @description
#' This function takes probability estimates from TransPhylo and runs a
#' linear regression on them, or logistic regression if a probability cutoff
#' to label an individual as an infection source is specified.
#' 
#' @param method A character string specifying the type of regression to run.
#' Can be "linear", "logistic", or "bayesian_logistic_misclass".
#' @param cleaned_data A data frame containing the covariates to use in the
#' regression. Must contain a column "SampleID" with the sample IDs.
#' @param prob_source A data frame of probabilities output by TransPhylo with columns
#' `SampleID` and `prob_source` (as output by `run_TransPhylo`), indicating
#' the probability of each sample being the infection source.
#' @param prob_cutoff A numeric specifying the probability cutoff to use for
#' logistic regression. If NULL, logistic regression will not be run.
#' @param sensitivity Numeric. Only used if `method` is "bayesian_logistic_misclass".
#' @param specificity Numeric. Only used if `method` is "bayesian_logistic_misclass".
#' @param ... Additional arguments to pass to Stan if `method` is
#' "bayesian_logistic_misclass". Iterations and chains default to 2000 and 4
#' respectively.
#'
#' @returns An `lm` object, a `glm` object, or a `stanfit`.
#' 
#' @importFrom rstan stan extract
#'
#' @export

regression <- function(method = c("logistic", "linear",
                                  "bayesian_logistic_misclass"),
                       cleaned_data, prob_source, prob_cutoff = NULL,
                       sensitivity = NULL, specificity = NULL,
                       ...) {
  # We will assume cleaned_data is a data frame with only the
  # covariates we want to use in the regression, and that the
  # sample ids are given in the column "SampleID"
  cleaned_data <- merge(cleaned_data, prob_source, by = "SampleID")
  cleaned_data <- cleaned_data[, !names(cleaned_data) %in% "SampleID"]
  # Drop rows with NA values
  cleaned_data <- cleaned_data[complete.cases(cleaned_data), ]
  if (method == "linear") {
    return(lm(prob_source ~ ., data = cleaned_data))
  }
  if (method == "logistic") {
    if (is.null(prob_cutoff)) {
      stop("Must provide a probability cutoff for logistic regression")
    }
    cleaned_data$tp_source <- as.factor(cleaned_data$prob_source > prob_cutoff)
    return(glm(tp_source ~ ., data = cleaned_data, family = "binomial"))
  }
  if (method == "bayesian_logistic_misclass") {
    if (is.null(prob_cutoff)) {
      stop("Must provide a probability cutoff for logistic regression")
    }
    if (is.null(sensitivity) || is.null(specificity)) {
      stop("Must provide sensitivity and specificity for Bayesian logistic regression")
    }
    # check that all columns are numeric
    if (!all(sapply(cleaned_data, is.numeric))) {
      stop("All columns in cleaned_data must be numeric for Bayesian logistic regression")
    }
    cleaned_data$tp_source <- as.numeric(cleaned_data$prob_source > prob_cutoff)
    covariates <- cleaned_data[, !(names(cleaned_data) %in% c("prob_source", "tp_source"))]
    # Run Isaac's Bayesian "homebrew": logistic regression with misclassification
    model_objects <- list(N = nrow(cleaned_data),
                          K = ncol(cleaned_data) - 2,
                          z = cleaned_data$tp_source,
                          X = covariates,
                          specificity = specificity,
                          sensitivity = sensitivity)
    stanfile <- system.file("Stan", "logistic_misclass.stan",
                            package = "beast2tpPipeline")
    args <- list(...)
    if (!"iter" %in% names(args)) {
      args$iter <- 2000
    }
    if (!"chains" %in% names(args)) {
      args$chains <- 4
    }
    args$file <- stanfile
    args$data <- model_objects
    stanfit <- do.call(rstan::stan, args)
    samples <- rstan::extract(stanfit)
    colnames(samples$beta) <- colnames(covariates)
    return(samples)
  }
}