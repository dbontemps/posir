
### utils ###

#' Put a path (directory) and a file name together in a string
#'
#' @param path_dir the path (string) of the directory.
#' @param file_without_path the file name (string, without path).
#'
#' @keywords internal
file_with_path_name <- function(path_dir, file_without_path) {
  fileN <- file_without_path
  if (path_dir != "") {
    if (substring(path_dir, nchar(path_dir)) == "/") {
      if (substr(file_without_path, 1, 1) == "/") fileN <- substring(file_without_path, 2)
    } else if (substr(file_without_path, 1, 1) != "/") fileN <- paste("/", file_without_path, sep = "")
    fileN <- paste(path_dir, fileN, sep = "")
  }
  return(fileN)
}

#' Gets the path to the "inst" sub-directory of posir
#'
#' @keywords internal
get_inst_posir_path <- function() {
  curwd <- getwd()
  if (substring(curwd, nchar(curwd)-5) == "/posir") {
    curwd <- paste(curwd, "inst", sep = "/")
    if (dir.exists(curwd)) return(curwd)
  }
  return(find.package("posir"))
}

#' Find the path to a file in the SavedOutputs sub-directory
#'
#' In the [posir] tree, the SavedOutputs sub-directory is used for testing;
#' it contains mainly .RData files and some .txt ones,
#' which each store the desired returns of some function,
#' or a toy dataset to be passed as argument.
#'
#' When not installed, SavedOutputs is a subdirectory of "posir/inst/";
#' when installed, SavedOutputs is a subdirectory of "posir/".
#'
#' find_Outputs_file() is a dev tool function,
#' probably destined to later disappear or at least be rewritten.
#' For now, it only works on my computers...
#'
#' @param OutputFN name (string, without path) of the saved file.
#' @param do_warn if FALSE, no warning is issued if the file is not found.
#'
#' @return the path to the saved file.
#'
#' @examples
#' \dontrun{
#' mypathQ <- find_Outputs_file("QTable_1040_10080.txt")
#' if (!file.exists(mypathQ)) log_n_stop("not found saved file QTable_1040_10080.txt")
#' }
#' @keywords internal
find_Outputs_file <- function(OutputFN, do_warn = TRUE) {
  res <- paste(get_inst_posir_path(), "SavedOutputs", sep = "/")
  if (OutputFN == "" || substr(OutputFN, 1, 1) == "/") {
    res <- paste(res, OutputFN, sep = "")
  } else {
    res <- paste(res, OutputFN, sep = "/")
  }
  if (do_warn && !file.exists(res)) {
    log_warn(paste("found no file", res))
    log_debug("in find_Outputs_file().")
  }
  return(res)
}

### Random simulations ###

#' Random generation following a centered or symmetrized Pareto distribution
#'
#' Simulate following a Pareto distribution, then symmetrize it or center it.
#'
#' If the symmetrized parameter is TRUE,
#' * takes n Pareto random variables;
#' * translates them so that their support is the nonnegative real numbers;
#' * then symmetrizes them (using a Bernoulli variable).
#'
#' Else, takes n Pareto random variables and centers them.
#'
#' @param n number of generated i.i.d. random variables.
#' @param location location parameter(s) of the Pareto distribution
#' (before centering or symetrizing), with default value 1.
#' @param shape shape parameter(s) of the Pareto distribution
#' (before centering or symetrizing), with default value 3.
#' @param symmetrized boolean(s) indicating whether the Pareto variables
#' are to be symmetrized or centered, with default value FALSE (centered).
#'
#' @return A vector of simulated centered or symmetrized Pareto random variables.
#' @importFrom EnvStats rpareto
#' @importFrom stats rbinom
#' @export
#'
#' @examples
#' rCenteredPareto(10, shape = c(1, 3), symmetrized = TRUE)
#' # Simulates 10 symmetrized Pareto,
#' # with location 1,
#' # and shape alternating between 1 and 3.
#'
rCenteredPareto <- function(n, location = 1, shape = 3, symmetrized = FALSE) {
  auxsym <- rep(symmetrized, length.out = n)
  auxsh <- rep(shape, length.out = n)
  auxloc <- rep(location, length.out = n)
  return(
    (rpareto(n, location, shape)
     - ifelse(auxsym, auxloc, auxsh * auxloc / (auxsh - 1))) # centering if not symmetrized, else just translating to R+ before symetrizing
    * ifelse(auxsym, 2 * rbinom(n, size = 1, prob = .5) - 1, 1) # symetrizing
  )
}

