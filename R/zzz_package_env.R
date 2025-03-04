# Package environment
pkg.env <- new.env(parent = emptyenv()) # pourquoi cet argument ? FIXME
# Comment passser des éléments d'environnement aux futures dans la parallélisation ?

## logging
pkg.env$cur.log.dir <- NULL
pkg.env$cur.log.file <- NULL
pkg.env$cur.do.tee <- FALSE

# POSIR processes quantiles tables
pkg.env$qtables <- list(NULL,NULL)

## Suite ne fonctionne pas FIXME
# Pour le moment :
#logger::log_threshold(logger::DEBUG, namespace = "posir")
#local_log_init()
# A terme :
#logger::log_threshold(logger::INFO, namespace = "posir")
#logger::log_appender(logger::appender_console, namespace = "posir")
