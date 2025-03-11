
#' Load the quantiles table of a POSIR process from a file
#'
#' The columns correspond to the delta parameter of the POSIR process.
#' Due to legacy code, the row correspond to a queue probability parameter
#' (alpha = 1-p). Both have to be sorted in decreasing order,
#' if not the function log a warning and returns NULL.
#'
#' @param NameF file name (to load from).
#' @return A table of Quantile Estimations for a POSIR process.
#'
#' @examples
#' \dontrun{
#' mypathQ <- paste(find.package("posir"),
#'   "/extdata/SavedOutputs/QTable_1040_10080.txt",
#'   sep = ""
#' )
#' read_quantiles(mypathQ)
#' }
#'
#' @keywords internal
read_quantiles <- function(NameF) {
  if (!file.exists(NameF)) {
    log_warn(
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

#' Load the pre-computed quantiles tables for 1D and 2D POSIR processes
#'
#' @keywords internal
load_quantiles <- function() {
  base_path <- paste(get_inst_posir_path(),
                     "/extdata/Table_quantiles_", sep="")
  for(i in 1:2) {
    pkg.env$qtables[[i]] <- read_quantiles(paste(
      base_path, toString(i), "D.txt", sep=""))
  }
}


#' Quantile distribution of POSIR processes
#'
#' Extrapolates a quantile from a previous table of (estimated) quantiles.
#'
#' @param alpha vector of queue probabilities;
#' usually a vector p of probabilities is used in such function,
#' but due to legacy code, here alpha == 1 - p.
#' @param delta vector of delta parameters;
#' when d==1, the delta parameters have to be in the interval \[0.005, 1\];
#' with d==2, the delta parameters have to be in the interval \[0.01, 1\];
#' in both cases, the estimated quantiles for delta very small are to be taken with caution.
#' @param d dimension: currently 1 and 2 are supported.
#' @param Qtable a table of (estimated) quantiles for a POSIR process;
#' if NULL, then the pre-computed tables supplied with the package are used;
#' otherwise, Qtable must be a matrix,
#' whose colnames indicate a grid of delta values,
#' whose rownames indicate a grid of probability values,
#' and where both are sorted in decreasing order.
#'
#' @returns A table of quantiles.
#'
#' @examples
#' \dontrun{
#' mypathQ <- paste(find.package("posir"),
#'   "/extdata/SavedOutputs/QTable_1040_10080.txt",
#'   sep = ""
#' )
#' posir_quantiles(NULL, .25, Qtable = read_quantiles(mypathQ))
#' }
#'
#' @keywords internal
posir_quantiles <- function(alpha, delta, d=1, Qtable=NULL) {
  if(is.null(Qtable)) {
    if(d != 1 && d!=2) log_n_stop("Currently unsupported dimension")
    if(is.null(pkg.env$qtables[[d]])) load_quantiles()
    Qtable <- pkg.env$qtables[[d]]
    if(is.null(Qtable))
      log_n_stop("The pre-compiled quantiles table could not be loaded")
  }
  alphagrid <- as.numeric(rownames(Qtable))
  if(is.null(alpha)) alpha <- alphagrid
  deltagrid <- as.numeric(colnames(Qtable))
  if(is.null(delta)) delta <- deltagrid
  # Without interpolation:
  #return(Qtable[sapply(p,toString), sapply(delta,toString)])
  # With linear interpolation:
  aux <- function(x,grid,gridname) {
    l <- length(grid)
    m <- length(x)
    p1 = p2 = 1:m
    c1 <- rep(1, m)
    c2 <- rep(0, m)
    for(i in 1:m) {
      p1[i] <- sum(grid>=x[i])
      p2[i] <- l+1-sum(grid<=x[i])
      if(p1[i]==0 || p2[i]>l)
        log_n_stop(paste("Out of range", gridname, "parameter in posir_quantiles()."))
      if(p2[i]!=p1[i]) {
        v1 <- grid[p1[i]]
        v2 <- grid[p2[i]]
        c1[i] <- (x[i]-v2)/(v1-v2)
        c2[i] <- (v1-x[i])/(v1-v2)
      }
    }
    return(list(p1=p1,p2=p2,c1=c1,c2=c2))
  }
  rowpos <- aux(alpha, alphagrid, "probability")
  colpos <- aux(delta, deltagrid, "delta")
  Q <- Qtable[rowpos$p1,colpos$p1]*rowpos$c1%*%t(colpos$c1) +
    Qtable[rowpos$p2,colpos$p1]*rowpos$c2%*%t(colpos$c1) +
    Qtable[rowpos$p1,colpos$p2]*rowpos$c1%*%t(colpos$c2) +
    Qtable[rowpos$p2,colpos$p2]*rowpos$c2%*%t(colpos$c2)
  colnames(Q) <- delta
  rownames(Q) <- alpha
  return(Q)
}

#' Quantile distribution of POSIR processes
#'
#' Extrapolates a quantile from the pre-computed estimated POSIR quantiles.
#'
#' @param p vector of probabilities.
#' @inheritParams posir_quantiles
#'
#' @returns A table of quantiles.
#' @export
#'
#' @examples
#' qposir(c(0.9,0.95,0.99), .25)
qposir <- function(p, delta, d=1) { posir_quantiles(1-p, delta, d) }


