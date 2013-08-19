\name{sosspp-package}
\alias{sosspp-package}
\alias{sosspp}
\docType{package}
\title{The SosSpp Package}
\description{
The package \pkg{sosspp} is a collection of statistical tools for the 
analysis of homogeneous and inhomogeneous spatial point patterns, 
on top of the \pkg{spatstat} package.
It accompanies the paper \emph{Inhomogeneous spatial point processes with hidden
second-order stationarity} by Ute Hahn and Eva B. Vedel Jensen.
}
\details{
The package \pkg{sosspp} deals with the statistical analysis 
of homogeneous and inhomogeneous point patterns that are supposed to originate 
from a reweighted, retransformed or rescaled second-order stationary point process. 
For a definition of these classes of spatial point processes, see the paper Hahn & 
Jensen (2013).

The \pkg{hidden2ndorder} package is build on top of the package \pkg{spatstat}, 
and makes use of its data formats. \pkg{sosspp} currently provides 
\itemize{
  \item estimates of the template \eqn{K}- and \eqn{L}-functions and of 
         the directional \eqn{K}-function
  \item model tests of the type of hidden second-order stationarity.
}
}
\section{Reweighted second-order stationary point processes}{
  Methods for this class of models are contained in the \pkg{spatstat}-package;
  the template \eqn{K}-function is estimated by \code{\link{Kinhom}}.
}
\section{Retransformed second-order stationary point processes}{
  To estimate the template \eqn{K}-function of a retransformed point process, 
  apply the estimate of the stationary \eqn{K}-function from \pkg{spatstat}, 
  \code{\link{Kest}}, to the backtransformed point pattern. General transformations
  of point patterns are obtained by the function
  \code{\link{coordTransform.ppp}}.
}
\section{Rescaled second-order stationary point processes}{
  The function \code{\link{Klocscaled}} replaces the \pkg{spatstat}-function 
  \code{\link{Kscaled}} for the estimation of the rescaled template \eqn{K}-function.
}
\author{
  Ute Hahn, \email{ute@imf.au.dk}
}
\references{
Hahn, U. and Jensen, E. B. V. (2013)
  Inhomogeneous spatial point processes with hidden second-order stationarity.
}
\keyword{ package }
\seealso{
  The \code{\link[spatstat:spatstat-package]{spatstat}} package
}
\examples{
# analysis of the scholtzia data set, same intensity estimate as in the article
data(scholtzia)
data(intensities)
# compare estimates of reweighted and rescaled template K-function
meanintens <- npoints(scholtzia) / area.owin(scholtzia)
r.ls <- seq(0,1.25,0.002)
r.inhom <- r.ls / sqrt(meanintens)
# obsolete! plot(Linhom(scholtzia, lambda = scholtzia.intens, renormalise = FALSE, 
#     r = r.inhom), main = "L.inhom Scholtzia")
# plot(Llocscaled(scholtzia, lambda = scholtzia.intens, renormalise = FALSE, 
#    r = r.ls), main = "L.locscaled Scholtzia")

# the testing
W1 <- owin(c(0, 22), c(0, 5.48), unitname=c("metre", "metres"))
W2 <- owin(c(0, 22), c(5.48, 22), unitname=c("metre", "metres"))
qwin1 <- quadrats(W1, nx = 4, ny = 1)
qwin2 <- quadrats(W2, nx = 2, ny = 2)

par(mfrow = c(2, 2))
plot(Kinhom(scholtzia[W1], lambda = scholtzia.intens, renormalise = FALSE, 
    r = r.inhom, correction="iso"), main = "K.inhom Scholtzia")
plot(Kinhom(scholtzia[W2], lambda = scholtzia.intens, renormalise = FALSE, 
    r = r.inhom, correction="iso"), add = TRUE)
# plot(Klocscaled(scholtzia[W1], lambda = scholtzia.intens, renormalise = FALSE, 
#     r = r.ls, correction="iso"), main = "K.locscaled Scholtzia")
# plot(Klocscaled(scholtzia[W2], lambda = scholtzia.intens, renormalise = FALSE, 
#     r = r.ls, correction="iso"), add = TRUE)


}
