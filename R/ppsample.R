# ppsubsamples


#'@title Subsample a point pattern
#@aliases ppsubsample
#'@description Subsample a point pattern on quadrats, generating a list of
#'point patterns (\code{\link{sostyppp}} objects)
#'@param pp object of class \code{\link{sostyppp}} or spatstat-class \code{\link[spatstat]{ppp}},
#'the point pattern to be subsampled
#'@param quads optional a list of objects of spatstat-class \code{\link[spatstat]{owin}} or a
#'spatstat - \code{\link[spatstat]{tess}} object
#'@param ... arguments passed to \code{spatstat::\link[spatstat]{quadrats}}
#'or an object of spatstat-class \code{\link[spatstat]{tess}}, the quadrats for subsampling \code{pp}.
#or a named list containing lists of \code{owin}s or \code{tess} objects, see the Details.
#'@return a \code{ppsample} object: a list of objects of same class as the input \code{pp}.
#'@details If \code{quads} is not given, the quadrats are determined by
#'\code{spatstat::\link[spatstat]{quadrats}}.
#'
#'The \code{ppsample} can be plotted, and number of points or estimated intensity
#'can be retrieved by methods \code{npoints} and \code{intensity}.
#'@export
#

ppsubsample <- function (pp, quads = NULL, ...)
{
  # names throw warnings in plyr, if used later.
  # Remove all names that are not list names from quads
  unnamelist <- function(l)
{
  if (!is.list(l)) names(l) <- NULL
  if (is.list(l)) for (i in seq_along(l)) l[[i]] <- unnamelist(l[[i]])
  return(l)
}

  if (is.null(quads)) quads <- quadrats(pp, ...)
  if("tess" %in% class(quads)) quads <- tiles(quads)
  nonamequads <- unnamelist(quads)
  ppsample <- lapply(nonamequads, function(w) pp[w])
#  npts <- sapply(ppsample, npoints)
#  areas <- sapply(ppsample, area.owin)
  attr(ppsample, "parentwindow") <- pp$window
#  attr(ppsample, "npoints") <- npts too dangerous if manipulated
#  attr(ppsample, "area") <- areas
#  attr(ppsample, "intensity") <- npts / areas
  if(!("ppsample" %in% class(ppsample))) class(ppsample) <- c("ppsample", class(ppsample))
  return(ppsample)
}

#' Check whether an object is a sample of point patterns
#'
#' Checks if an object belongs to class \code{"ppsample"}.
#'
#' @param x any \code{R} object
#' @return \code{TRUE} if \code{x} belongs to class \code{"ppsample"}, otherwise \code{FALSE}.
#' @export
# @seealso \code{\link{sostyppp.object}} for details on the class.
# @author Ute Hahn,  \email{ute@@imf.au.dk}

is.ppsample <- function(x) inherits(x, "ppsample")


#'Estimated intensities of a point pattern sample
#'
#'The estimated intensities for all point patterns contained in the data set.
#'
#'@param X a sample of point patterns, object of class \code{ppsample}
#'@param ... ignored, only for compatibility with generic method in \code{spatstat}
#'@return A numeric vector of estimated intensities.
#'@details
#'The intensities are estimated as the number of points divided by the area.
#'@S3method intensity ppsample
#'@method intensity ppsample
#@export
# @author Ute Hahn,  \email{ute@@imf.au.dk}
#'@seealso \code{\link{ppsubsample}} for creating \code{ppsample} objects,
#'\code{\link{npoints.ppsample}} for obtaining the number of points,
#'\code{\link[spatstat]{intensity}} for spatstats generic function.
#@examples

intensity.ppsample <- function(X, ...) sapply(X, npoints)/ sapply(X, area.owin)


#'Number of points in a point pattern sample
#'
#'Retrieves the number of points for all point patterns contained in the data set.
#'
#'@param x a sample of point patterns, object of class \code{ppsample}
#'@return integer vector, the number of points in each pattern contained in \code{X}
#'@S3method npoints ppsample
#'@method npoints ppsample
#@export
# @author Ute Hahn,  \email{ute@@imf.au.dk}
#'@seealso \code{\link{ppsubsample}} for creating \code{ppsample} objects,
#'\code{intensity.ppsample} for estimating the empirical intensity,
#'\code{\link[spatstat]{npoints}} for spatstat's generic function.
#@exampl

npoints.ppsample <- function(x) sapply(x, npoints)

#' Backtransformed a list of point pattern
#'
#' Backtransform all point patterns in a sample, as well as the parent window.
#'
#' @param X a sample of retransformed s.o.s. typed point patterns, object of class \code{{ppsample}}.
#' @return a \code{ppsample} object consisiting of homogeneous s.o.s. typed point patterns.
#' @details The parent window is also retransformed. For more details, see the function
#' \code{\link{backtransformed}} for single point patterns.
#' @S3method backtransformed ppsample
#' @method backtransformed ppsample
#' @export backtransformed.ppsample
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
# @examples
# bronzetra <- retransformed(bronzefilter, "gradx")
#'bronzesample <- ppsubsample(bronzetra, quadrats(bronzetra, nx=6, ny = 3))
#'plot(bronzesample, use.marks = FALSE)
#'plot(backtransformed(bronzesample), use.marks = FALSE)

backtransformed.ppsample <- function(X)
{
  stopifnot(is.ppsample(X))
  if (length(X) == 0) stop ("does not contain anything")
  X1 <- X[[1]]
  stopifnot(is.sostyppp(X1))
  stopifnot(hasType(X1, type= "t"))
  W <- attr(X, "parentwindow")
  if (!is.null(W)){
      sostinfo <- attr(X1, "sostinfo")
      W <- coordTransform(W,
                          trafoxy = sostinfo$backtransform,
                          invtrafoxy = sostinfo$transform,
                          subdivideBorder = !is.rectangle(W) | is.null(sostinfo$gradient))
  }
  Y <- lapply(X, backtransformed)
  attr(Y, "parentwindow") <- W
  if(!("ppsample" %in% class(Y))) class(Y) <- c("ppsample", class(Y))
  return(Y)
}



#' @title Extract subset of a point pattern sample
#' @aliases [.ppsample
#' @rdname Subset_ppsample
#' @description Subsample a point pattern sample, retaining information about
#' the original pattern's window.
#'
#' @return A \code{ppsample} object.
# @S3method [ ppsample
#' @method [ ppsample
#' @export
#' @param x a list of point patterns, object of class \code{"ppsample"}.
#' @param i subset index.
# @param j,drop ignored.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}


"[.ppsample" <- function(x, i) {
  y <- NextMethod()
  attr(y, "parentwindow") <- attr(x, "parentwindow")
  class(y) <- class(x)
  y
}