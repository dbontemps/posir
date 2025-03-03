
#' Estimate effective error levels from simulated trajectories and quantiles
#'
#' For each discretisation point \eqn{\delta} of the processes,
#' compute the proportion of trajectories that in \eqn{\delta}
#' are bigger than the corresponding quantile.
#'
#' Nothing forces us to only use this function with the 1D or 2D POSIR process!
#'
#' @param Zt transposed matrix of simulated trajectories (in columns),
#' for a discretized random process.
#' @param Q quantile table (matrix) of the limit random process.
#' @param precN how many decimals for the effective error level estimations
#' (used in particular in exporting the table to a file).
#' @seealso [extract_quantiles()], [run_simu()], [compute_error_levels()]
#'
#' @return A table (matrix) of estimated effective error levels.
#' @importFrom foreach foreach %do%
#' @export
#'
#' @examples
#' N <- 1000
#' dfs <- 3:5
#' alphagrid <- c(.5, .1, .05, .01)
#' m <- max(dfs)
#' Z <- matrix(stats::rnorm(N * m), N, m)^2
#' Zt <- as.matrix(apply(Z, 1, cumsum))[dfs, ] # apply transpose
#' rownames(Zt) <- dfs
#' Q <- as.matrix(cbind(sapply(dfs, function(df) {
#'   stats::qchisq(1 - alphagrid, df)
#' })))
#' colnames(Q) <- dfs
#' rownames(Q) <- alphagrid
#' extract_error_levels(Zt, Q, 3)
extract_error_levels <- function(Zt, Q, precN) {
  # Comparaison des quantiles
  j <- NULL
  levs <- foreach(j = 1:nrow(Q), .combine = "rbind", .multicombine = TRUE) %do% {
    round(rowMeans(Zt > Q[j, ]), precN)
  }
  colnames(levs) <- colnames(Q)
  rownames(levs) <- rownames(Q)
  return(levs)
}

#' Export effective error levels tables to files
#'
#' @param levs table of effective error levels.
#' @param NameF file name (to export to).
#' @param NameFQ file name of the quantile table used in estimation.
#' @param Ntraj number of trajectories
#' used in the Monte-Carlo estimation of the effective error levels.
#' @inheritParams simulationDir
#' @seealso [extract_error_levels()], [compute_error_levels()], [simulationDir()]
#'
#' @return NULL
#' @importFrom utils write.table
#'
#' @keywords internal
write_error_levels <- function(levs, NameF, NameFQ, Ndis, Ntraj, NameDis, d) {
  write("######\n#", NameF)
  write(
    paste("# Table of Effective Error Levels, \n",
          "# with alpha in rownames et delta in colnames, \n",
          "# for a ", toString(d), "D discretized POSIR process. \n",
          sep = ""
    ),
    NameF,
    append = TRUE
  )
  write(paste("#\n### Dimention", toString(d), "###\n#"), NameF, append = TRUE)
  write(paste("# Quantiles table:\n# \"", NameFQ, "\"", sep = ""),
        NameF,
        append = TRUE
  )
  write("#\n# Simulation parameters\n#", NameF, append = TRUE)
  write(paste("# n =", toString(Ndis)), NameF, append = TRUE)
  write(paste("# Ntraj =", toString(Ntraj)), NameF, append = TRUE)
  write(paste("# Distribution:", NameDis), NameF, append = TRUE)
  write("#\n######", NameF, append = TRUE)
  write.table(levs, NameF, append = TRUE)
  return(NULL)
}

### Main function for effective error levels ###

#' Computes and save effective error levels for the POSIR process
#'
#' Computes Monte-Carlo estimations
#' of the effective error levels of confidence intervals
#' for a discretized, possibly non-Gaussian, 1D or 2D POSIR process,
#' then save them in a file.
#'
#' @param NdisQ discretisation parameter previously used in computing the quantiles.
#' @param deltagrid grid for \eqn{\delta} previously used in computing the quantiles.
#' @param posdelta vector of indexes to define a sub-grid for \eqn{\delta}.
#' The resulting sub-grid (subset of the grid) deltagrid\[posdelta\]
#' is used in [run_simu()] to discretize the POSIR process trajectories.
#' @param gridname name of the sub-grid for \eqn{\delta},
#' to be used in [simulationDir()]
#' for the path definition of the error levels table.
#' @param longdgrid_name name of the initial grid for \eqn{\delta},
#' to be used in [simulationDir()] to retrieve the path of the quantile table.
#' @param NameFQ file name (with path) of the quantile table.
#' If NameFQ is NULL, the file name is computed anew using [simulationDir()]
#' from d, NdisQ, longdgrid_name, and filenamebasisQ.
#' @param filenamebasisQ file name for the quantile table (without the path).
#' @param filenamebasis file name for the error levels table (without the path).
#' @inheritParams n_traj_simu
#' @inheritParams run_simu
#' @inheritParams write_batch
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @inheritParams batch_filename
#' @inheritParams extract_error_levels
#' @seealso [run_simu()], [extract_error_levels()], [write_error_levels()], [compute_quantiles()]
#'
#' @return A table of (simultaneous) effective error levels,
#' only if the error levels file does not already exist,
#' else the function returns NULL.
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' logger::log_info("Running examples for compute_error_levels()",
#'                  namespace = "posir")
#' sim_dir <- local_init_testing(do_multi = TRUE,
#'                               envir=sys.frame(sys.nframe()))
#' mypathQ <- paste(find.package("posir"),
#'   "/SavedOutputs/QTable_1040_10080.txt",
#'   sep = ""
#' )
#' compute_error_levels(
#'   1000, 50, 1040, 40, seq(7, 2, -1) / 8, 1:6, "de .875 à .25",
#'   "6_div_8", 10^6, sim_dir, mypathQ
#' )
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
compute_error_levels <- function(Ntot, Batchsize, NdisQ, Ndis, deltagrid, posdelta,
                                 gridname, longdgrid_name, maxmatrixsize,
                                 sim_root_dir = "",
                                 NameFQ = NULL,
                                 filenamebasisQ = "Table_quantiles.txt",
                                 filenamebasis = "Table_niveaux_effectifs.txt",
                                 NameBaseBatch = "Batch",
                                 rdistrib = rCenteredPareto,
                                 NameDis = "Pareto3",
                                 precN = floor(log10(Ntot) / 2) + 2,
                                 d = 1, ErLev = .001) {
  local_log_init(sim_root_dir)
  if (is.null(NameFQ)) {
    log_info("Computing anew the path to the table file of the quantiles")
    log_debug("in compute_error_levels().")
    aux <- paste(simulationDir(sim_root_dir, NdisQ, longdgrid_name, "Gaussian", d),
      filenamebasisQ,
      sep = "/"
    )
  } else {
    aux <- NameFQ
  }
  # if (!file.exists(aux)) {
  #   log_n_stop(
  #     paste("Not found quantile table", aux),
  #     "in compute_error_levels()."
  #   )
  # }
  # Q <- as.matrix(read.table(aux))
  # titles <- paste("X", sapply(deltagrid, toString), sep = "")
  # l <- length(deltagrid)
  # if (dim(Q)[2] != l || sum(colnames(Q) == titles) < l) {
  #   log_n_stop(
  #     "not matching delta grid in quantile table",
  #     "in compute_error_levels()."
  #   )
  # }
  # colnames(Q) <- deltagrid
  Q <- read_quantiles(aux)
  if(is.null(Q) || dim(Q)[2] != length(deltagrid) ||
     sum(colnames(Q) == deltagrid) < length(deltagrid)) {
    log_n_stop("Invalid quantiles table in compute_error_levels().")
  }
  #
  NameF <- paste(simulationDir(sim_root_dir, Ndis, gridname, NameDis, d),
    filenamebasis,
    sep = "/"
  )
  if (!check_table_file(NameF, deltagrid[posdelta], nrow(Q))) {
    Z <- run_simu(Ntot, Batchsize, Ndis, deltagrid[posdelta], gridname,
      maxmatrixsize, sim_root_dir, NameBaseBatch, rdistrib, NameDis,
      is_standard = FALSE, d, ErLev
    )
    log_debug("compute_error_levels(): now to the error levels!")
    Ntraj <- nrow(Z)
    levs <- extract_error_levels(t(Z), Q[, posdelta], precN)
    write_error_levels(levs, NameF, NameFQ, Ndis, Ntraj, NameDis, d)
    log_debug("compute_error_levels(): done.")
    return(levs)
  }
  return(NULL)
}
