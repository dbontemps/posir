

### Wrappers for logger functions ###

#' Wrap and hide [logger::log_debug()]
#'
#' @param message whatever you want to log.
#' @importFrom logger log_debug
#'
#' @keywords internal
log_debug <- function(message) {
  logger::log_debug(message, namespace = "posir")
}

#' Wrap and hide [logger::log_info()]
#'
#' @inheritParams log_debug
#' @importFrom logger log_info
#'
#' @keywords internal
log_info <- function(message) {
  logger::log_info(message, namespace = "posir")
}

#' Wrap and hide [logger::log_warn()]
#'
#' @inheritParams log_debug
#' @importFrom logger log_warn
#'
#' @keywords internal
log_warn <- function(message) {
  logger::log_warn(message, namespace = "posir")
}

#' Wrap and hide [logger::log_error()]
#'
#' @inheritParams log_debug
#' @importFrom logger log_error
#'
#' @keywords internal
log_error <- function(message) {
  logger::log_error(message, namespace = "posir")
}

### Other log functions ###

#' Log and raise an exception (stop) from R
#'
#' A message is logged at level ERROR, then stop(message) is called.
#' An optional additional message (details) can be given,
#' to be logged at level DEBUG.
#'
#' @param message text (string) to be logged and passed to stop().
#' @param m_details (string) optional additional details to be log at DEBUG level.
#' @family log and stop from R or C++
#'
#' @keywords internal
log_n_stop <- function(message, m_details = "") {
  log_error(message)
  # if(m_details=="") stop(message)
  # else {
  #   log_debug(m_details)
  #   stop(paste(mdetails, message, sep=": "))
  # }
  if (m_details != "") log_debug(m_details)
  stop(message)
}

### Setting up the logging system ###

#' Logging inside a POSIR simulation directory
#'
#' @param sim_dir the POSIR simulation directory
#' in which a log file is to be filled.
#' @param logFN the name (without path) of the log file.
#' @param do_tee if TRUE, logging is also done in the console.
#'
#' @importFrom logger log_appender appender_console appender_tee appender_file
#' @return the path to the log file.
#'
#' @examples
#' \dontrun{
#' logger::log_info("Running examples for log_in_dir()",
#'                  namespace = "posir")
#' fileName = "mylog.log"
#' sim_dir = local_init_testing(dirNameBase = "posir_test", logFN = fileName,
#'                               do_log=FALSE, envir=sys.frame(sys.nframe()))
#' olddir = pkg.env$cur.log.dir
#' oldfile = pkg.env$cur.log.file
#' olddotee = pkg.env$cur.do.tee
#' myfile = log_in_dir(sim_dir, fileName, TRUE)
#' logger::log_info(
#'   paste("We are logging in dir ", sim_dir,
#'     "/, in file ", fileName, "\n",
#'     sep = ""
#'   ),
#'   namespace = "posir"
#' )
#' # checking
#' cat(paste("contents of file", myfile, "\n"))
#' for (l in readLines(myfile)) cat(paste(l, "\n"))
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
#' log_in_dir(olddir, oldfile, olddotee)
#' }
#' @keywords internal
log_in_dir <- function(sim_dir, logFN, do_tee) {
  myfile = NULL
  if(is.null(sim_dir) || is.null(logFN)) {
    log_appender(appender_console, namespace = "posir")
  } else {
    myfile = file_with_path_name(sim_dir, logFN)
    if (sim_dir != "") {
      if (!dir.exists(sim_dir)) {
        log_debug(paste("log_in_dir(): creating directory.", sim_dir))
        dir.create(sim_dir)
      }
    }
    if (!file.exists(myfile)) file.create(myfile)
    log_debug(paste("Have set up log file", myfile, "in log_in_dir()."))
    if (do_tee) {
      log_appender(force(appender_tee(file = myfile)), namespace = "posir")
    } else {
      log_appender(force(appender_file(myfile)), namespace = "posir")
    }
  }
  pkg.env$cur.log.dir = sim_dir
  pkg.env$cur.log.file = logFN
  pkg.env$cur.do.tee = do_tee
  if(is.null(myfile)) {
    log_debug("Have set up logging in console")
    log_debug("in log_in_dir().")
  } else {
    log_debug("Initializing the new log system.")
  }
  return(myfile)
}

#' Initiate logging if it has not already been done
#'
#' local_log_init() sets up logging in the simulation directory,
#' [withr][withr::withr]-style.
#'
#' @inheritParams log_in_dir
#' @param envir  [environment] the environment to use for scoping.
#' @importFrom withr defer
#'
#' @keywords internal
local_log_init <- function(sim_dir=getwd(), logFN = "posir.log",
                           do_tee = TRUE, envir = parent.frame()) {
  olddir = pkg.env$cur.log.dir
  oldfile = pkg.env$cur.log.file
  olddotee = pkg.env$cur.do.tee
  if(sim_dir!=olddir || logFN!=oldfile || do_tee!=olddotee) {
    if(!is.null(olddir)) log_debug(paste("Old logging dir:",olddir))
    log_debug(paste("Setting up logging in file",
                    file_with_path_name(sim_dir, logFN),
                    "in local_log_init()."))
    log_in_dir(sim_dir, logFN, do_tee)
    defer({
      log_debug("Restoring previous logging environment.")
      log_in_dir(olddir, oldfile, olddotee)
    }, envir=envir)
  }
}
