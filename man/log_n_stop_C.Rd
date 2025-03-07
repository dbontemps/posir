% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{log_n_stop_C}
\alias{log_n_stop_C}
\title{Log and raise an exception (stop) from C++}
\usage{
log_n_stop_C(message, m_details)
}
\arguments{
\item{message}{text (string) to be logged and passed to stop().}

\item{m_details}{(string) optional additional details to be log at DEBUG level.}
}
\description{
A message is logged at level ERROR, then stop(message) is called.
An optional additional message (details) can be given,
to be logged at level DEBUG.
}
\examples{
\dontrun{
  tryCatch(log_n_stop_C("essai de log_n_top_C", "via tryCatch()"),
           error= function(e){log_info(paste("erreur :",e))},
           finally=function(e){log_info(paste("finally :",e))})
}
}
\seealso{
Other log and stop from R or C++: 
\code{\link{log_n_stop}()}
}
\concept{log and stop from R or C++}
\keyword{internal}
