#' @title Use Effective Sample Size to Check Mixing
#'
#' @description
#' For a BEAST2 run or a TransPhylo run, pull the effective sample size
#' for estimated parameters. Check that they are above some threshold.
#'
#' @param program Character string. Either "BEAST2" or "TransPhylo".
#' @param path_to_mcmc_log Character string. Path to the MCMC log file.
#' @param burn_in_fraction Numeric. Fraction of the MCMC chain to discard as
#' burn-in.
#' @param sample_interval Numeric. Interval at which samples were taken.
#' @param min_ess Numeric. Minimum effective sample size to consider acceptable.
#'
#' @returns NULL
#'
#' @import coda
#' @import tracerer
#' @import TransPhylo
#'
#' @export
ess_checks <- function(program = c("BEAST2", "TransPhylo"),
                       # for BEAST2 this is log file,
                       # for TP it's the resTransPhylo object
                       path_to_mcmc_log,
                       burn_in_fraction = 0.1,
                       sample_interval = 5000,
                       min_ess) {
  if (program == "BEAST2") {
    beast_log_full <- tracerer::parse_beast_tracelog_file(path_to_mcmc_log)
    beast_log <- tracerer::remove_burn_ins(beast_log_full,
                                           burn_in_fraction = burn_in_fraction)
    ess <- tracerer::calc_esses(beast_log, sample_interval = sample_interval)
  } else if (program == "TransPhylo") {
    tp_res <- readRDS(path_to_mcmc_log)
    ess <- lapply(tp_res, function(res) {
      class(res) <- "resTransPhylo"
      res <- TransPhylo::convertToCoda(res)
      coda::effectiveSize(res)
    })
    ess <- unlist(ess)
  }
  if (min(ess) < min_ess) {
    too_low <- ess[which(ess < min_ess)]
    nm <- basename(path_to_mcmc_log)
    message(paste0(nm, ":\n", "ESS is too low for ",
                   length(too_low), " parameters"))
    message(paste0("\t", "ESS: ", too_low))
  } else {
    message("All ESS are above minimum threshold")
  }
}