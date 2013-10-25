#' @useDynLib sostatpp
#'
# wishes: name statspatpp statspp sostatspat SOStatSpatPP spat2ostat spatsostat sostatpp sstatpp
# wishes: package options for plotting (lightcol=0.66)
#'
#'
#'@import fdnonpar
#'@import plutils
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
#'It implements methods from the papers \emph{"A studentized permutation test
#'for the comparison of spatial point patterns"} and \emph{"Inhomogeneous spatial
#'point processes with hidden second-order stationarity"} by Ute Hahn and
#'Eva B. Vedel Jensen.
#'
#'Hahn & Jensen (2013) introduce the term "hidden second-order stationarity"
#'as an umbrella for several classes of inhomogeneous models, viz.
#'reweighted, retransformed or rescaled second-order stationary point processes.
#'These model classes all come with their own sensible definition of an appropriate
#'inhomogenous version of second-order summary statistics, such as the \eqn{K}-function.
#'H&J (2013) use the notion of a "template version" for these second-order statistics.
#'
#'The \pkg{sostatpp} package extends the package \pkg{spatstat}, and
#'makes use of its data formats. \pkg{sostatpp} currently provides
#'\itemize{
#'  \item estimates of the template \eqn{K}- and \eqn{L}-functions and of
#'   the \eqn{\Delta K_{dir}}-function, described in H&J(2013)
#' \item a model test of the type of hidden second-order stationarity,
#' \item a test to compare the \eqn{K}-function of two point patterns, as in H(2012).
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
#'  \cr\code{\link{homogeneous}} \tab for plain, homogeneous s.o. stationarity.
#'  }
#'  A \code{sostyppp}-object can contain information on several types of s.o.
#'  stationarity. To retrieve the type(s) of an \code{sostyppp}, apply
#'  \tabular{ll}{
#'  \code{\link{hasType}} \tab to check if the object contains information for a particular type
#'  \cr\code{\link{currentType}} \tab to retrieve the s.o.s.-type that will be used for analysis
#'  }
#'@section Data manipulation and subsampling:
#'\subsection{Transformation}{
#' Arbitrary coordinate transformation on \pkg{spatstat} objects can be done with
#'     \code{\link{coordTransform}}, currently supporting
#'\tabular{ll}{
#'  \code{\link{coordTransform.im}} \tab transformation of \code{\link{im.object}s},
#'  \cr\code{\link{coordTransform.owin}} \tab transformation of \code{\link{owin.object}s},
#'  \cr\code{\link{coordTransform.ppp}} \tab transformation of \code{\link{ppp.object}s}.
#'}
#'
#'\code{\link{backtransformed}}  returnes the backtransform of a retransformed s.o.s. point pattern
#'or of a subsample of such a pattern
#'}
#'\subsection{Subsampling}{
#'\tabular{ll}{
#'\code{ppsubsample}\tab generates a quadrat subsample of a point pattern (object of type \code{ppsample}) 
#'\cr\code{plot.ppsample}\tab can be used to plot a subsample
#'\cr\code{twoquadsets}\tab to divide a set of quadrats into two sets that suitable 
#'for testing second-order stationarity
#'\cr\code{quadratsplot}\tab a plot to visualize two quadrat subsamples of a point pattern.
#'}
#'}
#'@section Statistical analysis:
#'\subsection{Second order summary functions}{
#'The following functions return objects of \pkg{spatstat}-class \code{fv}:
#' \tabular{ll}{
#' \code{\link{K.est}} \tab estimates the template \eqn{K}-function, according to type
#' \cr\code{\link{L.est}} \tab estimates the template \eqn{L}-function,
#' \cr\code{\link{DeltaKdir.est}} \tab estimates the \eqn{\Delta K_{dir}}-function
#' }
#' }
#' \subsection{Tests}{
#' \tabular{ll}{
#' \code{\link{sos.test}} \tab test type of s.o.s or local anisotropy,
#'    using \eqn{K}-function on subsamples
#' \cr\code{\link{twosample.K.test}} \tab comparison of the \eqn{K}-function
#'   on two samples of point patterns.
# \cr\code{\link{Kaniso.test}} \tab test of local anisotropy, using \eqn{\Delta K_{dir}}-function
#'}
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
#'  \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#'
#'
#'@keywords  package
#'@seealso  The \code{\link[spatstat]{spatstat}} package
#'@examples
#'## The data analysis section in H&J(2013) is covered by the following demos: ##
#'demo(bronzefilter)
#'demo(beilschmiedia)
#'demo(scholtzia)
NA
