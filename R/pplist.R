# ppsubsamples


#'@title Subsample a point pattern
#'@description Subsample a point pattern on quadrats, generating a list of 
#'point patterns (\code{\link{sostyppp}} objects)
#'@param pp object of class \code{\link{sostyppp}} or spatstat-class \code{\link{ppp}}, 
#'the point pattern to be subsampled
#'@param quads optional a list of objects of spatstat-class \code{\link{owin}} or a
#'spatstat - \code{\link{tess}} object
#'@param ... arguments passed to spatstat-function \code{\link{quadrats}}
#'or an object of spatstat-class \code{\link{tess}}, the quadrats for subsampling \code{pp}. 
#or a named list containing lists of \code{owin}s or \code{tess} objects, see the Details.
#'@return a \code{\link{ppsample}} object: a list of objects of same class as the input \code{pp}.
#'@details If \code{quads} is not given, the quadrats are determined by spatstat-function 
#'\code{\link{quadrats}}.
#'
#'The \code{ppsample} can be plotted, and number of points or estimated intensity
#'can be retrieved by methods \code{npoints} and \code{intensity}.
#'@export
#

quadratsubsample <- function (pp, quads = NULL, ...)
{
  # names throw warnings in plyr, if used later. 
  # Remove all names that are not list names from quads  
  unnamelist <- function(l)
{
  if (!is.list(l)) names(l) <- NULL
  if (is.list(l)) for (i in 1:length(l)) l[[i]] <- unnamelist(l[[i]])
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
  class(ppsample) <- c("ppsample", class(ppsample))
  return(ppsample)
}


#'Estimated intensities of a point pattern sample
#'
#'The estimated intensities for all point patterns contained in the data set.
#'
#'@param X a sample of point patterns, object of class \code{ppsample}
#'@return A numeric vector of estimated intensities.
#'@details
#'The intensities are estimated as the number of points divided by the area.
#'@S3method intensity ppsample
#'@method intensity ppsample 
#@export
# @author Ute Hahn,  \email{ute@@imf.au.dk}
#'@seealso \code{\link{ppsubsample}} for creating \code{ppsample} objects, 
#'\code{\link{npoints.ppsample}} for obtaining the number of points,
#'\code{\link{intensity}} for spatstats generic function.
#@examples

intensity.ppsample <- function(X) sapply(X, npoints)/ sapply(X, area.owin)


#'Number of points in a point pattern sample
#'
#'Retrieves the number of points for all point patterns contained in the data set.
#'
#'@param X a sample of point patterns, object of class \code{ppsample}
#'@return integer vector, the number of points in each pattern contained in \code{X}
#'@S3method npoints ppsample
#'@method npoints ppsample 
#@export
# @author Ute Hahn,  \email{ute@@imf.au.dk}
#'@seealso \code{\link{ppsubsample}} for creating \code{ppsample} objects, 
#'\code{intensity.ppsample} for estimating the empirical intensity, 
#'\code{\link{npoints}} for spatstats generic function.
#@examples

npoints.ppsample <- function(X) sapply(X, npoints)

