
### Selection of a temp directory ###

#' Select and create a temporary directory for simulations
#'
#' @param dirNameBase the name of a subdirectory in which run the simulations;
#' dirNameBase must not be "" nor ".", nor contain "..";
#' if subdirectories of subdirectories are used, only the last one may not exist.
#' @param create_if_needed if TRUE and if the planned directory does not exist,
#' it is created.
#' @return A path (string) in which simulations can be run.
#' @importFrom stringr str_detect
#'
#' @examples
#' \dontrun{
#' logger::log_info("Running examples for temp_root_dir()",
#'                  namespace = "posir")
#' sim_dir <- local_init_testing(dirNameBase = "posir_test",
#'                               envir=sys.frame(sys.nframe()))
#' myfile <- paste(sim_dir, "bidon.txt", sep = "/")
#' write("0", myfile)
#' print(file.exists(myfile))
#' # Last checking
#' if (dir.exists(sim_dir)) {
#'   log_error("something went wrong in temp_root_dir()")
#' }
#' withr::deffered_run()
#' }
#' @keywords internal
temp_root_dir <- function(dirNameBase = "POSIR_simu", create_if_needed=TRUE) {
  res <- tempdir()
  if (dirNameBase == "" || dirNameBase == ".") {
    log_n_stop("must use a sub-directory.", "in temp_root_dir().")
  }
  if (str_detect(dirNameBase, "\\.\\.")) {
    log_n_stop("no going back in the path!", "in temp_root_dir().")
  }
  res <- paste(res, dirNameBase, sep = "/")
  if (create_if_needed && !dir.exists(res)) dir.create(res)
  return(res)
}


### initializing temporary environment for POSIR simulation and tests ###

#' Initialize the directory for POSIR simulations
#'
#' local_init_testing() sets up the simulation directory
#' for testing various posir functions, [withr][withr::withr]-style.
#'
#' If do_temp is TRUE, a temporary directory is used
#' as parent of the simulation directory;
#' if do_log is TRUE, [local_log_init()] is called;
#' if do_multisession is TRUE, a future multisession plan is registered.
#'
#' @param do_temp if TRUE, a temporary directory is used.
#' @param do_log if TRUE, logging is set up in the simulation directory.
#' @param do_multi if TRUE, future::plan(multisession) is called.
#' @param do_seed if TRUE, a fixed (arbitrary) seed is set up.
#' @inheritParams temp_root_dir
#' @inheritParams log_in_dir
#' @inheritParams local_log_init
#' @importFrom future plan multisession
#' @importFrom withr defer local_rng_version local_seed
#' @seealso [withr::defer()], [withr::local_seed()].
#'
#' @return the path of the simulation directory for testing.
#' @export
#'
#' @examples
#' my_const_random = function() { local_init_testing(); rnorm(1) }
#' rnorm(1)
#' my_const_random()
#' rnorm(1)
#' my_const_random()
#' rnorm(1)
local_init_testing = function(do_temp = TRUE, dirNameBase = "POSIR_simu",
                              do_log = FALSE, logFN = "posir.log", do_tee = TRUE,
                              do_multi = FALSE, do_seed = TRUE,
                              envir=parent.frame()) {
  #envir = sys.frame(sys.nframe())
  # logging is last one, so that in the rest of the current logging system.
  if (do_seed) {
    local_rng_version("4.1.2", .local_envir = envir)
    local_seed(34781, .local_envir = envir, .rng_kind = "Mersenne-Twister",
               .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  }
  if (do_temp) {
    mypath <- temp_root_dir(dirNameBase, FALSE)
    if(dir.exists(mypath))
    {
      log_debug(paste("Will not delete already present testing directory",
                      mypath))
      log_debug("in local_init_testing().")
    }
    else {
      dir.create(mypath)
      defer(unlink(mypath, recursive = TRUE), envir = envir)
    }
  } else {
    mypath <- dirNameBase
    if (!dir.exists(dirNameBase)) {
      dir.create(dirNameBase)
      log_debug(paste("Will delete previously absent testing directory",
                      dirNameBase))
      log_debug("in local_init_testing().")
      defer(unlink(dirNameBase, recursive = TRUE), envir = envir)
    }
  }
  if (do_log) local_log_init(mypath, logFN, do_tee, envir = envir)
  if (do_multi) {
    prev_plan <- plan(multisession)
    defer(plan(prev_plan), envir = envir)
  }
  return(mypath)
}
