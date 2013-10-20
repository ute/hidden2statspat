#' Estimate summary function on quadrats
#'
#' Estimate  a  summary function on subsamples of a point pattern.
#'
#' @param pp point pattern, object of  class \code{\link{sostyppp}} or 
#' {spatstat}-class \code{\link{ppp}}
#' @param quads quadrats for subsampling \code{pp}. A \code{list} of objects
#' of spatstat-class \code{\link{owin}} or an object of spatstat-class \code{\link{tess}}.
#' @param fun the summary function to be applied
#' @param ... further arguments passed to \code{"fun"} 
#' @return An object of class \code{fvlist}, which is a list of 
#'   spatstat-\code{\link{fv}} objects.
# @seealso \code{\link{estK}}, \code{\link{estDeltaKdir}}
#    for estimation of the template \eqn{K}-function or of Delta K_dir.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @keywords robust
#' @keywords spatial
#' @keywords nonparametric
#' @keywords ts
#'@examples
#'pp <- rpoispp(100, win=square(4))
#'Ks <- fv.on.quadrats (pp, quadrats(pp, nx = 2, ny = 2))
#'lapply(Ks, plot)
fv.on.quadrats <- function(pp, quads, fun = estK, ...) {
  ppsample <- ppsubsample(pp, quads)
  fvl <- lapply(ppsample, fun, ...)
  class(fvl) <- c("fvlist", "list")
  fvl
}

#'@title Coerce a list of \code{fv} objects into a \code{fdsample}
#'@description Transform a list of spatstat function tables (\code{\link{fv}} objects) 
#'into an \code{\link{fdsample}} object, extracting one of the columns in each 
#'function table.
#'@param fvl list of spatstat function tables (\code{fv} objects),
#'@param valname character, the name of the column to be extracted. If not given, the
#'recommended value stored in attribute \code{"valu"} of the \code{fv} objects is 
#'used.
#'@return An \code{\link{fdsample}} object.
#'@details An \code{fdsample} stores function values in a matrix. The function
#'values in every column of the matrix are supposed to refer to the same argument
#'values. \code{as.fdsample.fvlist} does not interpolate values in \code{fvl}, 
#'therefore, the function tables in \code{fvl} must have identical argument values.
#'
#'The plot options \code{xlab} and \code{ylab} for the resulting \code{fdsample}
#'are retrieved from the entries in \code{fvl}.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @keywords robust
#' @keywords spatial
#' @keywords nonparametric
#' @keywords ts
#'@examples 
#'Ks <- replicate(4, list(Kest(rpoispp(100), r = seq(0, 0.2, 0.005))))
#'Ksample <- extract.fdsample(Ks, "iso")
#'plot(Ksample)
extract.fdsample  <- function(fvl, valname = NULL) {
  OK <- inherits(fvl, "fvlist") || (is.list(fvl) && all(sapply(fvl, is.fv)))
  if (!OK) 
    stop("not a list of fv")
  if (length(fvl) < 1)
    return(fdsample(numeric(0), numeric(0)))
  # spatstat's collapse.fv anyway needs arguments...
  if (!is.null(valname)){
    OK <- all(sapply(fvl, function(x) valname %in% attr(x, "dotnames")))
  }
  else {
    valname <- attr(fvl[[1]], "valu")
    OK <- all(sapply(fvl, function(x) valname %in% attr(x, "dotnames")))
  }
  if (!OK)
    stop("correction", corr, "not found in all members of fv-list")
  argu <- attr(fvl[[1]], "argu")
  args <- as.list(fvl[[1]])[[argu]]
  if (!all(sapply(fvl, function(x) as.list(x)[[argu]]==args)))
    stop("not all members of fv-list have the same argument values")
  fvals <- sapply(fvl, function(x) as.list(x)[[valname]])
  ylab <- attr(fvl[[1]], "ylab")
  if (is.language(ylab)) ylab <- as.expression(ylab)
  fdsample(args, fvals, xlab = argu, ylab = ylab)
}