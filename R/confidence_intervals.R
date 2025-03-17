
#' Check whether the arguments are within specification for 1D POSIR
#'
#' @param x a vector of observations to which a segmentation is applied.
#' @param breakpoints a vector of breakpoints positions
#' (the last row index before a breakpoint)
#' which must be given in increasing order.
#' @seealso [CI1Dparamvalidation()]
#'
#' @returns a boolean
#'
#' @keywords internal
CI1Ddatavalidation <- function(x, breakpoints) {
  if(is.numeric(x) && is.vector(x) &&
     (is.null(breakpoints) ||
      (is.numeric(breakpoints) && is.vector(breakpoints)))
     ) {
    n=length(x)
    l=length(breakpoints)+1
    aux_bkp=c(0, breakpoints, n)
    if(sum(sapply(1:l, function(y) {return(aux_bkp[y]==round(aux_bkp[y]) &&
                                           aux_bkp[y]<aux_bkp[y+1])})
    )==l) { return(TRUE) }}
  return(FALSE)
}

#' Check whether the arguments are within specification for 1D POSIR
#'
#' @param delta a lower bound on the relative length of the regions,
#' under which the confidence intervals are not given. Currently in dimension 1,
#' delta has to be in the interval \[0.005, 1\].
#' @param alpha the simultaneous error level for the confidence intervals.
#' @seealso [CI1Ddatavalidation()]
#'
#' @returns a boolean
#'
#' @keywords internal
CI1Dparamvalidation <- function(delta, alpha) {
  return(is.numeric(delta) && is.vector(delta) &&
           length(delta)==1 && delta<=1 && delta>=0.005 &&
           is.numeric(alpha) && is.vector(alpha) &&
           length(alpha)==1 && alpha<=0.5 && alpha>=0.001)
}

#' Check the validity of the elements of an object of class posir1DCI.
#'
#' @param y an object of class posir1DCI.
#' @seealso [new_posir1DCI()]
#'
#' @returns a boolean.
#' @export
#'
#' @examples
#' x=posir1DCI(rnorm(100),c(),0,-.196,.196)
#' validate_posir1DCI(x)
validate_posir1DCI <- function(y) {
  if(is.list(y) && CI1Ddatavalidation(y$data,y$breakpoints)) {
    l=length(y$breakpoints)+1
    if(is.numeric(y$lbs) && length(y$lbs)==l &&
       is.numeric(y$ubs) && length(y$ubs)==l &&
       (is.null(y$means) || (is.numeric(y$means) && length(y$means)==l))
       ) {
      for(i in 1:l) {
        if(is.na(y$lbs[i])) {
          if(!is.na(y$ubs[i])) return(FALSE)
        } else {
          if(is.na(y$ubs[i])) return(FALSE)
          if(y$ubs[i]<y$lbs[i]) return(FALSE)
        }
      }
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Create an object of class posir1DCI.
#'
#' @inheritParams CI1Ddatavalidation
#' @param means a vector of mean values
#' for each region obtained by the segmentation.
#' @param lbs a vector of lower bounds for the confidence intervals.
#' @param ubs a vector of upper bounds for the confidence intervals.
#'
#' @returns an object of class posir1DCI, which is a list with 4 fields:
#' * $data: a vector of data upon which a segmentation is applied;
#' * $breakpoints: a vector of breakpoints
#'   (the last row index before a breakpoint);
#' * $means: the means in each intervals (if given, otherwise set to NULL).
#' * $lbs a vector for the lower bounds
#'   of the confidence intervals in each segment;
#' * $ubs a vector for the upper bounds
#'   of the confidence intervals in each segment.
#' @export
#'
#' @examples
#' posir1DCI(rnorm(100),c(),0,-.196,.196)
posir1DCI <- function(x=0, breakpoints=numeric(), means=NULL,
                      lbs=0, ubs=0) {
  res=list(data=x, breakpoints=breakpoints, means=means, lbs=lbs, ubs=ubs)
  class(res)="posir1DCI"
  if(!validate_posir1DCI(res)) log_n_stop("Invalid arguments in new_posir1DCI().")
  return(res)
}

#' Create a new object of class posir1DCI.
#'
#' @export
#' @seealso [posir1DCI()]
#'
new_posir1DCI <- function() { return(posir1DCI(0,numeric(),NULL,0,0)) }


#' Plot the 1D confidence intervals obtained by POSIR method.
#'
#' @param x a posir1DCI object.
#' @inheritParams base::plot
#' @seealso [posir1DCI()]
#'
#' @export
#' @importFrom graphics abline lines
#'
#' @examples
#' x=rnorm(50)+c(rep(1,15), rep(2,35))
#' bkp=c(10,20)
#' y=confidence_intervals_1D(x,bkp,.1,.005)
#' plot(y)
plot.posir1DCI <- function(x,...) {
  if(!validate_posir1DCI(x))
    log_n_stop("Invalid arguments in plot.posir1DCI().")
  n=length(x$data)
  l=length(x$breakpoints)+1
  aux_bkp=c(0, x$breakpoints, n)
  plotmeans=!is.null(x$means)
  plot(1:n,x$data,main="Simultaneous confidence intervals for the means",
       xlab="", ylab="",
       col = "#33333333", pch = 19, cex = 0.3,
       ...)
  for(i in 1:(l-1)) abline(v=x$breakpoints[i]+.5, col="red")
  for(i in 1:l) {
    xx=c(aux_bkp[i]+.5,aux_bkp[i+1]+.5)
    if(plotmeans) lines(xx, rep(x$means[i],2), col="red")
    if(!is.na(x$lbs[i])) {
      lines(xx, rep(x$lbs[i],2),
            col="red", lty=2)
      lines(xx, rep(x$ubs[i],2),
            col="red", lty=2)
    }
  }
}

#' Compute the simultaneous confidence intervals
#' for a segmentation applied to a numeric vector
#' using 1D POSIR method.
#'
#' @inheritParams CI1Ddatavalidation
#' @inheritParams CI1Dparamvalidation
#' @param homoscedatic a boolean which indicates whether
#' the estimation of the standard deviation is common across the segments
#' (when true), or is done separately in each segment (when false).
#'
#' @returns an object of class posir1DCI, which contains
#' simultaneous confidence intervals for the mean value of data x
#' in each region delimited by breakpoints.
#'
#' @export
#'
#' @examples
#' x=rnorm(50)+c(rep(1,15), rep(2,35))
#' bkp=c(10,20)
#' confidence_intervals_1D(x,bkp,.1,.005)
confidence_intervals_1D <- function(x, breakpoints, delta, alpha,
                                    homoscedatic = TRUE) {
  if(!CI1Ddatavalidation(x, breakpoints) ||
     !CI1Dparamvalidation(delta, alpha))
    log_n_stop("invalid argument in confidence_intervals_1D().")
  n=length(x)
  l=length(breakpoints)+1
  if(l>=n) log_n_stop("too many breakpoints in confidence_intervals_1D().")
  aux_bkp=c(0, breakpoints, n)
  someNA=FALSE
  lmin=max(2,ceiling(n*delta)) #
  q=posir_quantiles(alpha, delta, d=1)
  lbs=ubs=rep(NA,l)
  moy=lens=s=rep(0,l)
  for(i in 1:l) {
    cur = x[(aux_bkp[i]+1):aux_bkp[i+1]]
    lens[i]=length(cur)
    moy[i] = mean(cur)
    s[i] = sd(cur)
  }
  if(homoscedatic) {
    aux=c()
    for(i in 1:l) aux=c(aux,rep(moy[i],lens[i]))
    s=rep(sqrt(sum((x-aux)^2)/(n-l)),l)
  }
  for(i in 1:l) {
    if(lens[i]>=lmin) {
      r=q*s[i]/sqrt(lens[i])
      lbs[i]=moy[i]-r
      ubs[i]=moy[i]+r
    }
    else someNA=TRUE
  }
  if(someNA) log_warn("In confidence_intervals_1D(): some intervals in the segmentation were not long enough, so the corresponding confidence intervals could not be computed.")
  return(posir1DCI(x, breakpoints, moy, lbs, ubs))
}

