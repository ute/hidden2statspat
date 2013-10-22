#' @useDynLib sostatpp
#'
# wishes: name statspatpp statspp sostatspat SOStatSpatPP spat2ostat spatsostat sostatpp sstatpp
# wishes: package options for plotting (lightcol=0.66)
#'
#'
#'@import fdnonpar
#'@import plottools
#'@import plyr
#'@import spatstat
#'
#'
#'@name sostatpp-package
#'@aliases sostatpp-package sostatpp
#'@docType package
#'@title The Sostatpp Package
#'@description
#'The package \pkg{sostatpp} is a collection of statistical tools for the
#'analysis of homogeneous and inhomogeneous spatial point patterns, with
#'particular focus on (hidden) second-order stationarity. It supplements the
#'very comprehensive \pkg{spatstat} package.
#'@details
#'The package \pkg{sostatpp} deals with the statistical analysis of homogeneous
#'and inhomogeneous second-order stationary point patterns.
#'It implements methods from the papers \emph{A studentized permutation test
#'for the comparison of spatial point patterns} and \emph{Inhomogeneous spatial
#'point processes with hidden second-order stationarity} by Ute Hahn and
#'Eva B. Vedel Jensen.
#'
#'Hahn & Jensen (2013) introduce the term "hidden second-order stationarity"
#'as an umbrella of several classes of inhomogeneous models, viz.
#'reweighted, retransformed or rescaled second-order stationary point processes.
#'These model classes all come with their own sensible definition of an inhomogenous
#'"template" version to second-order summary statistics, such as the \eqn{K}-function.
#'
#'The \pkg{sostatpp} package extends the package \pkg{spatstat}, and
#'makes use of its data formats. \pkg{sostatpp} currently provides
#'\itemize{
#'  \item estimates of the template \eqn{K}- and \eqn{L}-functions and of
#'   the \eqn{\Delta K_{dir}}-function, described in H&J(2013)
#' \item model tests of the type of hidden second-order stationarity that can
#' also be used to compare two point patterns, as in H(2012).
#' }
#'
#' The three different types of second-order stationarity show many similarities.
#' In order to stress this aspect, the statistical analysis is unified by assigning
#' the stationarity type to the data beforehand. To this end, \pkg{sostatpp}
#' declares a class \code{\link{sostyppp}} of second-order stationarity typed point
#' patterns.
#'
#'
#' @section Assigning and retrieving the type of second-order stationarity:
#'  Any spatstat or sostatpp point pattern (object of class \code{ppp} or \code{sostyppp})
#'  can be converted into a second-order stationarity typed point pattern
#'  (object of class \code{sostyppp}), using
#'  \tabular{ll}{
#'  \code{\link{reweighted}} \tab for reweighted s.o.stationarity
#'  \cr\code{\link{retransformed}} \tab for retransformed s.o. stationarity
#'  \cr\code{\link{rescaled}} \tab for locally rescaled s.o.stationarity
#'  \cr\code{\link{ashomogeneous}} \tab for plain, homogeneous s.o. stationarity.
#'  }
#'  A \code{sostyppp}-object can contain information on several types of s.o.
#'  stationarity. To retrieve the type(s) of an \code{sostyppp}, apply
#'  \tabular{ll}{
#'  \code{\link{has.type}} \tab check for a particular type
#'  \cr\code{\link{currenttype}} \tab the s.o.s.-type that will be used for analysis
#'  }
#'@section Statistical analysis:
#'\subsection{Second order summary functions}{
#' \tabular{ll}{
#' \code{\link{K.est}} \tab estimates the \eqn{K}-function, according to type
#' \cr\code{\link{L.est}} \tab estimates the \eqn{L}-function,
#' \cr\code{\link{DeltaKdir.est}} \tab estimates the \eqn{\Delta K_{dir}}-function
#' }
#' \subsection{Tests}{
#' \tabular{ll}{
#' \code{\link{Kpermute.test}} \tab test type of s.o.s or local anisotropy,
#'    using \eqn{K}-function on subsamples
# \cr\code{\link{Kaniso.test}} \tab test of local anisotropy, using \eqn{\Delta K_{dir}}-function
#'}
#'}
#'}
#' @section Coordinate transformation:
#'   Arbitrary coordinate transformation on \pkg{spatstat} objects can be done with
#'     \code{\link{coordTransform}}, currently supporting
#'\tabular{ll}{
#'  \code{\link{coordTransform.im}} \tab transformation of \code{\link{im.object}s},
#'  \cr\code{\link{coordTransform.owin}} \tab transformation of \code{\link{owin.object}s},
#'  \cr\code{\link{coordTransform.ppp}} \tab transformation of \code{\link{ppp.object}s}.
#'}
#'\tabular{ll}{
#'\code{\link{backtransformed}} \tab returnes the backtransform of a retransformed s.o.s. point pattern
#'\cr
#'    }
#'@author Ute Hahn, \email{ute@@imf.au.dk}
#'
#'@references
#'Hahn, U. (2012) A studentized permutation test for the comparison of spatial point  patterns.
#'\emph{Journal of the American Statistical Association} \strong{107} (498), 754--764.
#'
#'Hahn, U. and Jensen, E. B. V. (2013)
#'  Inhomogeneous spatial point processes with hidden second-order stationarity.
#'  \emph{CSGB preprint} 2013-7.
#'  \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#'
#'
#'@keywords  package
#'@seealso  The \code{\link[spatstat]{spatstat}} package
#'@examples
#'#### analysis of the scholtzia data set ################################
#'
#' data(scholtzia)
#' data(intensities) # same intensity estimate as in the article
#' # compare estimates of reweighted and rescaled template K-function
#'  r.ls <- seq(0,1.25,0.002)
#'  r.inhom <- r.ls / sqrt(intensity(scholtzia))
#'
#'############### Exploring and testing hypotheses #############################
#'
#'######################### subsample ####################################
#'
#'medy <- median(scholtzia$y)
#'W1 <- owin(c(0, 22), c(0, medy), unitname=c("metre", "metres"))
#'W2 <-  owin(c(0, 22), c(medy, 22), unitname=c("metre", "metres"))
#'quads1 <- tiles(quadrats(W1, nx=4, ny=1))
#'quads2 <- tiles(quadrats(W2, nx=2, ny=2))
#'
#'sscholtzia <- rescaled(scholtzia, intensity = scholtzia.intens)
#'
#'ssample <- list(hi = ppsubsample(sscholtzia, quads1),
#'                lo = ppsubsample(sscholtzia, quads2))
################ estimate K on subsamples ###############################
#
## Ks1 <- estOnQuadrats(scholtzia, fun = K.est, type="s", lambda=scholtzia.intens,
##                      quads = quads1, rmax = 1.25)
#
##Ks2 <- estOnQuadrats(scholtzia, fun = K.est, type="s", lambda=scholtzia.intens,
##                      quads = quads2, rmax = 1.25)
#
######## plot the estimated K_functions ##################################
#
#plot(Ks1, col = "red", light = .5,
#     ylim=c(0, 6),  main=" rescaled, original intensity")
#plot(Ks2, col = "blue", light = .5, add=TRUE)
#
##########  test s.-o. rescaled stationarity ###########################
#
#Kpermute.test(scholtzia, type="s", quads1 = quads1, quads2 = quads2,
#                  rmax = 1.25, lambda=scholtzia.intens, use.tbar = TRUE)
NA
