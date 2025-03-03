# Package environment
pkg.env <- new.env()

## logging
pkg.env$cur.log.dir <- NULL
pkg.env$cur.log.file <- NULL
pkg.env$cur.do.tee <- FALSE

# POSIR processes quantiles tables
pkg.qtables <- list(NULL,NULL)
load_quantiles()

## Suite ne fonctionne pas FIXME
# Pour le moment :
#logger::log_threshold(logger::DEBUG, namespace = "posir")
#local_log_init()
# A terme :
#logger::log_threshold(logger::INFO, namespace = "posir")
#logger::log_appender(logger::appender_console, namespace = "posir")
