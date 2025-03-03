
#' Load the quantiles table of a POSIR process from a file
#'
#' Both the delta parameter (columns) and the probaility parameter (rows)
#' have to be sorted in decreasing order,
#' if not the function log a warning and returns NULL.
#'
#' @param NameF file name (to load from).
#' @return A table of Quantile Estimations for a POSIR process.
#'
#' @keywords internal
read_quantiles <- function(NameF) {
  if (!file.exists(NameF)) {
    log_info(
      paste("Not found quantile table", NameF,
      "in read_quantiles().")
    )
    return(NULL)
  }
  Q <- as.matrix(read.table(NameF))
  deltagrid <- as.numeric(substring(colnames(Q),2))
  colnames(Q) <- deltagrid
  alphagrid <- as.numeric(rownames(Q))
  aux <- function(grid,gridname) {
    l <- length(grid)
    if(l==0 || grid[1]>1 || grid[l]<=0 || sum(grid[-l]<=grid[-1])>0){
      log_warn(paste("Invalid grid for", gridname, "in read_quantiles()."))
      return(TRUE)
    }
    return(FALSE)
  }
  if(aux(deltagrid,"delta") || aux(alphagrid, "alpha")) return(NULL)
  return(Q)
}


