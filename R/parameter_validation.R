
#' Validate a grid against a discretisation parameter
#'
#' Check whether a vector is a decreasing subsequence
#' with values in \eqn{{ k/Ndis : 1 \le k \le Ndis }}.
#'
#' An exception is raised if the grid is invalid.
#' Else, the vector of integers deltagrid*Ndis is returned.
#'
#' @param Ndis discretisation parameter (of type integer).
#' @param deltagrid the grid to be checked (of type vector).
#' @param ErLev tolerance level against which a number is considered an integer.
#'
#' @return The vector deltagrid*Ndis, only if the grid is validated.
#'
#' @examples
#' \dontrun{
#' check_grid(30, seq(10, 2, -1) / 10, 10^(-3)) # valid grid
#' }
#' \dontrun{
#' check_grid(5, seq(10, 2, -1) / 10, .01) # invalid grid
#' }
#' @keywords internal
check_grid <- function(Ndis, deltagrid, ErLev) {
  aux <- deltagrid * Ndis
  pos_grid <- round(aux)
  l <- length(deltagrid)
  if (sum(abs(aux - pos_grid) < ErLev) < l) {
    log_n_stop(
      "delta grid and discretization are not compatible",
      "in check_grid()."
    )
  }
  if (sum(deltagrid > 0) < l || sum(deltagrid <= 1) < l ||
      sum(pos_grid[2:l] >= pos_grid[1:(l - 1)]) > 0) {
    log_n_stop("invalid delta grid", "in check_grid().")
  }
  return(pos_grid)
}


#' Validate previously saved tables (matrices) in files
#'
#' Check whether a file contains a table (matrix) of assigned dimensions,
#' and with the correct column names.
#'
#' If the file exists but does not contain a valid table, it is deleted.
#'
#' @param NameF filename (with path, of type string).
#' @param deltagrid grid for \eqn{\delta}.
#' @param rowsnb rows number.
#' @seealso [colnames()]
#'
#' @return A boolean, with value TRUE if a valid table has been found in the file.
#' @importFrom utils read.table
#'
#' @examples
#' \dontrun{
#' logger::log_info("Running examples for check_table_file()",
#'                  namespace = "posir")
#' sim_dir <- local_init_testing(dirNameBase = "posir_test",
#'                               envir=sys.frame(sys.nframe()))
#' newone <- paste(sim_dir, "mytemptable.txt", sep = "/")
#' X <- matrix(1, 50, 6)
#' colnames(X) <- seq(7, 2, -1) / 8
#' write.table(X, newone)
#' check_table_file(newone, seq(7, 2, -1) / 8, 50)
#' check_table_file(newone, seq(7, 2, -1) / 8, 100)
#' check_table_file(newone, seq(7, 2, -1) / 8, 50)
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
#' }
#' @keywords internal
check_table_file <- function(NameF, deltagrid, rowsnb) {
  if (file.exists(NameF)) {
    l <- length(deltagrid)
    Q <- read.table(NameF)
    titles <- paste("X", sapply(deltagrid, toString), sep = "")
    if (sum(dim(Q) == c(rowsnb, l)) < 2 || sum(colnames(Q) == titles) < l) {
      log_info(paste("Removing invalid table file \"", NameF, "\".", sep = ""))
      file.remove(NameF)
      return(FALSE)
    } else {
      log_debug(paste("Validated table file \"", NameF, "\".", sep = ""))
      return(TRUE)
    }
  } else {
    log_debug(paste("Found no table file \"", NameF, "\".", sep = ""))
    return(FALSE)
  }
}
