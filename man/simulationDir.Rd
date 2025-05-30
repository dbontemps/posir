% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batches.R
\name{simulationDir}
\alias{simulationDir}
\title{Define the directory for a POSIR simulation}
\usage{
simulationDir(sim_root_dir, Ndis, gridname, NameDis, d)
}
\arguments{
\item{sim_root_dir}{root directory of the simulations. Both "" and "."
stand for the current working directory.}

\item{Ndis}{discretisation parameter of the POSIR process.}

\item{gridname}{name of the grid for \eqn{\delta}.}

\item{NameDis}{name, or short description (of type string),
of the distribution of the i.i.d. random variables
that are summed in the POSIR process.}

\item{d}{dimension parameter of the POSIR process (1 or 2).}
}
\value{
Something similar to
"1D_Discretization_1080_Gaussian_grid_de .875 à .25.txt".
}
\description{
Compute the relative path (directory name) for stocking the batch files
in a POSIR simulation.
}
\details{
If the directory does not exist, it is created.
}
\examples{
\dontrun{
sim_dir <- local_init_testing(dirNameBase = "posir_test",
                              envir=sys.frame(sys.nframe()))
simulationDir(sim_dir, 1080, "de .875 à .25", "Gaussian", 1)
# some code
withr::deferred_run(envir=sys.frame(sys.nframe()))
}
}
\seealso{
\code{\link[=batch_filename]{batch_filename()}}
}
\keyword{internal}
