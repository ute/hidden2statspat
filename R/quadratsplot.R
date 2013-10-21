#'@title Plot a point pattern with test quadrats
#'@description A plot of a point pattern, with subdivision into quadrats. 
#'Several sets of quadrats can be plotted in the same plot, distinguished by 
#'different background colour.
#'@param x point pattern
#'@param quads the quadrat sets, a \emph{named} list of  lists of windows (\code{"\link{owin}"}) 
#' or tessellations (\pkg{spatstat}-object \code{"\link{tess}"})
#'@param styles optional a \emph{named} list of \code{\link{style}} lists, see Details,
#'@param ... further plot parameters
#'@details
#'Background colour and coulour intensity of the quadrats are determined by the
#'parameters \code{col.win} and \code{alpha.win}, see \code{\link{plot.sostyppp}}.
#'These plot parameters can be given as \code{\link{style}} lists or explicitely.
#'Both styles and explicite parameters can be given as named lists, the names are 
#'then matched with the names of the quadrat sets.
#'Quadrats that have not been assigned
#'@export
#'@author Ute Hahn,  \email{ute@@imf.au.dk}
#'@examples
#'# quadrat plot for the beilschmiedia pattern
#'beiquads <- quadshilo(bei, nx = 8, ny = 4, minpoints = 30)
#'beistyle <- list(hi = style(col.win = "red"), lo = style(col.win = "blue"))
#'
#'quadratsplot(bei, beiquads, beistyle, alpha.win = 0.2, main = "Beilschmiedia: quadrats")
#'
#'# does the same colours:
#'
#'quadratsplot(bei, beiquads, col.win = list(hi = "red", lo = "blue"),
#'  alpha.win = 0.2, main = "Beilschmiedia: quadrats", pch = 16, cex = .3)

quadratsplot <- function(x, qsets, styles = NULL, ...){
  x <- as.sostyppp(x, type = "none") # to be on the safe side
  ppsamples <- mapply(ppsubsample, quads = qsets, MoreArgs = list(pp = x))
  dotargs <- style(...)
  splot(x, matching(dotargs, plot.ppp, .graphparams))
  lplot(ppsamples, styles, ..., allinone = TRUE, add = TRUE)
}

