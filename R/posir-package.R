#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' posir
#'
#' @name posir-package
# #' @docType package
#' @import rlang
#' @importFrom glue glue
#' @importFrom lifecycle deprecated
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @useDynLib posir, .registration = TRUE
#'
#' @details
#' # Presentation
#' With posir, you can
#' * Simulate trajectories of 1D or 2D (discretized) POSIR processes,
#' * Estimate the quantiles under the standard Gaussian distribution,
#' * Estimate the effective error levels of confidence intervals,
#'   for discretized, possibly non-Gaussian, POSIR processes.
#'
#' # Available (exported) functions
#' * [compute_error_levels()]
#' * [compute_quantiles()]
#' * [extract_error_levels()]
#' * [extract_quantiles()]
#' * [n_traj_simu()], et ses versions \code{R}
#'   [n_traj_simu_1D_C()] et [n_traj_simu_2D_C()]
#' * [rCenteredPareto()]
#' * [run_simu()]
#'
#' # POSIR processes
#' ## Definition
#' See paper \insertCite{BaBoNe24}{posir}.
#'
#' ## Usage of POSIR processes
#' The POSIR process is used in a POSI-like scheme to obtain
#' asymptotically valid, post region selection, simultaneous confidence intervals
#' for the means of a signal over a set of regions,
#' in the presence of a white noise.
#'
#' # Simulations results
#' The folder "extdata/scripts/" in the posir package installation path
#' contains the codes used to perform simulations
#' for the paper \insertCite{BaBoNe24}{posir},
#' while the file "extdata/results.zip" contains the obtained results.
#'
#' @references
#'   \insertAllCited{}
#'
## usethis namespace: end
NULL
