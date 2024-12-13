#' @title Run Linear or Logistic Regression with TransPhylo Results
#'
#' @description
#' This function takes probability estimates from TransPhylo and runs a
#' linear regression on them, or logistic regression if a probability cutoff
#' to label an individual as an infection source is specified.
#' 
#' @param method A character string specifying the type of regression to run.
#' Can be "linear" or "logistic".
#' @param cleaned_data A data frame containing the covariates to use in the
#' regression. Must contain a column "SampleID" with the sample IDs.
#' @param prob_source A vector of probabilities output by TransPhylo, indicating
#' the probability of each sample being the infection source.
#' @param prob_cutoff A numeric specifying the probability cutoff to use for
#' logistic regression. If NULL, logistic regression will not be run.
#'
#' @returns An `lm` or `glm` object.
#'
#' @export

regression <- function(method = c("logistic", "linear"), cleaned_data,
                       prob_source, prob_cutoff = NULL) {
  # We will assume cleaned_data is a data frame with only the
  # covariates we want to use in the regression, and that the
  # sample ids are given in the column "SampleID"
  cleaned_data <- merge(cleaned_data, prob_source, by = "SampleID")
  cleaned_data <- cleaned_data[, !names(cleaned_data) %in% "SampleID"]
  if (method == "linear") {
    return(lm(prob_source ~ ., data = cleaned_data))
  }
  if (method == "logistic") {
    if (is.null(prob_cutoff)) {
      stop("Must provide a probability cutoff for logistic regression")
    }
    cleaned_data$prob_source <- as.factor(cleaned_data$prob_source > prob_cutoff)
    return(glm(prob_source ~ ., data = cleaned_data, family = "binomial"))
  }
}