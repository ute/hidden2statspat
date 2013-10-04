#' @useDynLib sostatpp
#' 
# wishes: name statspatpp statspp sostatspat SOStatSpatPP spat2ostat spatsostat sostatpp sstatpp
# wishes: package options for plotting (lightcol=0.66)
#'
#'
#'@import fdnonpar
# require(fdnonpar)
#'
#'
#'@name sostatpp-package
#'@aliases sostatpp-package sostatpp
#'@docType package
#'@title The Sostatpp Package
#'@description
#'The package \pkg{sostatpp} is a collection of statistical tools for the 
#'analysis of homogeneous and inhomogeneous second-order stationary spatial 
#'point patterns, on top of the \pkg{spatstat} package.
#'@details
#'The package \pkg{sostatpp} deals with the statistical analysis of homogeneous 
#'and inhomogeneous point patterns that are supposed to originate from a 
#'reweighted, retransformed or rescaled second-order stationary point process. 
#'For a definition of these classes of spatial point processes, see the paper Hahn 
#'& Jensen (2013).
#'It implements methods from the papers \emph{A studentized permutation test 
#'for the comparison of spatial point patterns} and \emph{Inhomogeneous spatial 
#'point processes with hidden second-order stationarity} by Ute Hahn and 
#'Eva B. Vedel Jensen.
#'
#'The \pkg{sostatpp} package is build on top of the package \pkg{spatstat}, and 
#'makes use of its data formats. \pkg{sostatpp} currently provides 
#'\itemize{
#'  \item estimates of the template \eqn{K}- and \eqn{L}-functions and of 
#'   the directional \eqn{K}-function
#' \item model tests of the type of hidden second-order stationarity that can 
#' also be used to compare two point patterns.
#' }
#' 
#' The three different types of second-order stationarity show many similarities.
#' In order to stress this aspect, the statistical analysis is unified by assigning
#' the stationarity type to the data beforehand. To this end, \pkg{sostatpp} 
#' declares a class \code{sostpp} of second-order stationarity typed point 
#' patterns, see \link{sostpp.object}.
#' 
#' 
#' @section {Assigning and retrieving the type of second-order stationarity}{
#'  Any spatstat or sostatpp point pattern (object of class \code{ppp} or \code{sostpp}) 
#'  can be converted into a second-order stationarity typed point pattern 
#'  (object of class \code{sostpp}), using 
#'  \tabular{ll}{
#'  \code{\link{reweighted}} \tab for reweighted s.o.stationarity 
#'  \cr\code{\link{retransformed}} \tab for retransformed s.o. stationarity 
#'  \cr\code{\link{rescaled}} \tab for locally rescaled s.o.stationarity 
#'  \cr\code{\link{ashomogeneous}} \tab for plain, homogeneous s.o. stationarity.
#'  }
#'  A \code{sostpp}-object can contain information on several types of s.o. 
#'  stationarity. To retrieve the type(s) of an \code{sostpp}, apply
#'  \tabular{ll}{
#'  \code{\link{has.type}} \tab check for a particular type 
#'  \cr \code{\link{currenttype}} \tab the s.o.s.-type that will be used for analysis
#'  }
#'}
#'@section {Statistical analysis}{
#'\subsection{Second order summary functions}{
#' \tabular{ll}{    
#' \code{\link{Khidden}} \tab estimates the \eqn{K}-function, according to type
#' \cr\code{\link{Lhidden}} \tab estimates the \eqn{L}-function,
#' \cr\code{\link{DeltaKDir}} \tab estimates the \eqn{\Delta K_{dir}}-function
#' }
#' }
#' \subsection{Tests}{
#' \tabular{ll}{
#' \code{\link{Kpermute.test}} \tab test type of s.o.s, using \eqn{K}-function on subsamples
#' \cr\code{\link{Kaniso.test}} \tab test of local anisotropy, using \eqn{\Delta K_{dir}}-function
#'}
#'}
#'}
#' @section{Coordinate transformation}{
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
#'}
#'@author Ute Hahn, \email{ute@@imf.au.dk}
#'
#'@references 
#'Hahn, U. (2012) A studentized permutation test for the comparison of spatial point  patterns.
#'\emph{Journal of the American Statistical Association} \strong{107} (498), 754--764.
#'
#'Hahn, U. and Jensen, E. B. V. (2013)
#'  Inhomogeneous spatial point processes with hidden second-order stationarity.
#'  \emph{CSGB preprint} 2013-7. 
#'
#'
#'@keywords  package 
#'@seealso  The \code{\link[spatstat:spatstat-package]{spatstat}} package
#'@examples 
#'#### analysis of the scholtzia data set ################################
#'
#' data(scholtzia)
#' data(intensities) # same intensity estimate as in the article
#' # compare estimates of reweighted and rescaled template K-function
#'  meanintens <- npoints(scholtzia) / area.owin(scholtzia)
#'  r.ls <- seq(0,1.25,0.002)
#'  r.inhom <- r.ls / sqrt(meanintens)
#'   
#'############### Testing hypotheses #######################################
#'
#'medy <- median(scholtzia$y)
#'W1 <- owin(c(0, 22), c(0, medy), unitname=c("metre", "metres"))
#'W2 <-  owin(c(0, 22), c(medy, 22), unitname=c("metre", "metres"))
#'quads1 <- tiles(quadrats(W1, nx=4, ny=1))
#'quads2 <- tiles(quadrats(W2, nx=2, ny=2))
#'
#'Ks1 <- estOnQuadrats(scholtzia, fun = estK, type="s", lambda=scholtzia.intens, 
#'                      quads = quads1, rmax = 1.25)
#'                      
#'Ks2 <- estOnQuadrats(scholtzia, fun = estK, type="s", lambda=scholtzia.intens, 
#'                      quads = quads2, rmax = 1.25)
#'plot(Ks1, col = "red", ylim=c(0, 6),  main=" rescaled, original intensity")
#'plot(Ks2, col = "blue", add=T, )
NA
