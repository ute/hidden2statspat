#'@title Plot a point pattern with test quadrats
#'@description A plot of a point pattern, with subdivision into quadrats.
#'Several sets of quadrats can be plotted in the same plot, distinguished by
#'different background colour.
#'@param x point pattern
#'@param qsets the quadrat sets, a \emph{named} list of  lists of windows (\code{"\link{owin}"})
#' or tessellations (\pkg{spatstat}-object \code{"\link[spatstat]{tess}"})
#'@param styles optional a \emph{named} list of \code{\link{simplist}}s, see Details,
#'@param ... further plot parameters
#'@param main title to be put above the plot. Default: tries to guess title from call.
#'@param backtransformed logical, if \code{TRUE} backtransform  pattern and quadrats.
#'Requires pattern of sos-type \code{"t"}.
#'#'@details
#'Background colour and coulour intensity of the quadrats are determined by the
#'parameters \code{col.win} and \code{alpha.win}, see \code{\link{plot.sostyppp}}.
#'These plot parameters can be given as \code{\link{simplist}}s or explicitely.
#'Both styles and explicite parameters can be given as named lists, the names are
#'then matched with the names of the quadrat sets.
#'Quadrats that have not been assigned
#'@export
#'@author Ute Hahn,  \email{ute@@imf.au.dk}
#'@examples
#'# quadrat plot for the beilschmiedia pattern
#'beiquads <- twoquadsets(bei, nx = 8, ny = 4, minpoints = 30)
#'beistyle <- list(hi = simplist(col.win = "red"), lo = simplist(col.win = "blue"))
#'
#'quadratsplot(bei, beiquads, beistyle, alpha.win = 0.2, main = "Beilschmiedia: quadrats")
#'
#'# does the same colours:
#'
#'quadratsplot(bei, beiquads, col.win = list(hi = "red", lo = "blue"),
#'  alpha.win = 0.2, main = "Beilschmiedia: quadrats", pch = 16, cex = .3)

quadratsplot <- function(x, qsets, styles = NULL, ..., main = NULL,
                        backtransformed = FALSE) {
  # qsets should be a list
  stopifnot(is.list(qsets), is.ppp(x))

  if (is.null(main)) {
    main <- deparse(substitute(x))
    if (backtransformed) main <- paste(main, ", backtransformed", sep = "")
  }
  if(!is.sostyppp(x)) x <- as.sostyppp(x, type = "none")
  ppsamples <- lapply(qsets, function(sa) ppsubsample(x, sa))

  if (backtransformed && identical(currentType(x), "t")){
    x <- backtransformed(x)
    ppsamples <- lapply(ppsamples, backtransformed.ppsample)
  }


  # plot original pattern, just to have unused points also in the plot
  do.call("plot.ppp", c(list(x),
        matching(simplist(..., main = main), plot.ppp, .graphparams)))
  # don't throw warnings because of extra arguments in the styles; 
  okstyles <- lapply(styles, matching, plot.ppp, plot.sostyppp, .graphparams)
  # force lty = "solid"
  lplot(ppsamples, okstyles, ..., lty = "solid", allinone = TRUE, add = TRUE)
}