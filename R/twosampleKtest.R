#'@title Two sample K-permutation-test
#'@description Perform a permutation test to compare K-functions estimated on
#'two point pattern samples
#' Returns an object of class \code{Ktest}
#'that also can be plotted, see the details.
#'@param x,y the point pattern samples to be compared, objects of class \code{ppsample}
#'@param r0 numeric, the upper integration limit, see Details,
#'@param rlen optional, number of steps for numerical integration, defaults to 256,
#see Details,
#'@param Kfun optional \code{function}, the \eqn{K}-function to be used,
#'  either \code{\link{estK}} (default) or \code{\link{estDeltaKdir}},
#'@param ... optional parameters for the function \code{Kfun}. To speed up 
#'calculations, it is recommended to give an explicite \code{correction} argument.
#'@param use.tbar logical, if true, a modified test statistic is used, see Details,
#'@param nperm number of random permutations, see Details. Defaults to 1000.
#@param noTest optional logical, if \code{TRUE}, no test is run, only the point 
#pattern subsamples and \eqn{K}-function estimates are returned.
#'
#'@details The function tests if the \eqn{K}-functions estimated on the
#'point pattern samples have the same mean. 
#'\subsection{What the test does, and details on the parameters}{
#'The \eqn{K}-function, or \eqn{\Delta K_{dir}}, is estimated on all patterns in
#'the two samples, and the resulting two samples of estimated K-functions are 
#'compared by a permutation test. 
#'The test statistic is the integral over a squared Welch-t-statistic,
#'\deqn{T=\int_0^{r_max}\frac{(K_1(r)-K_2(r))^2}{s_1^2(r)/m_1 + 
#'   s_2^2(r)/m_2} d r}{T = integral [ (K_1(r)-K_2(r))^2 / (s_1^2(r)/m_1 + 
#'   s_2^2(r)/m_2)],} 
#'where \eqn{K_1(r)} and \eqn{K_2(r)} are the group means, and
#'\eqn{s_1^2(r), s_2^2(r)} are within group variances at a fixed argument \eqn{r}.
#'The integral spans an interval \eqn{[0, r_{max}]}. It is approximated by the mean
#'of the integrand over all \code{rlen} values of \eqn{r}, multiplied by the length
#'of the integral, i.e., by \code{r0}.
#'
#'A variant of the test statistic, \eqn{\bar T}{Tbar}, uses the variance stabilized
#'function \eqn{K(r)/r} instead of \eqn{K(r)} and replaces the denominator in
#'the integrand with  \eqn{mean (s_1^2(r)/m_1 + s_2^2(r)/m_2)}. To use this variant
#'instead of the original statistc \eqn{T}, let \code{use.tbar = TRUE}
#'
#'The \emph{p}-value is obtained by permutation across the groups; the number of
#'permutations is specified by \code{nperm}. If \code{nperm = NULL}, 
#'the exact test with all permutations is used (combinations, for symmetry reasons). 
#'This may cause memory or computing time issues. 
#'If \code{nperm} is given as an integer, the permutations are sampled randomly, 
#'unless \code{nperm} is larger than the number of disjoint combinations. 
#'In that case, the exact version is applied, see \code{\link{tL2.permtest}}.
#'
#'To test against differences in local anisotropy,  H\&J (2013) propose to use 
#'the \eqn{\Delta K_{dir}}-function 
#'instead of the isotropic \eqn{K}-function. 
#'For this variant of the test, let \code{Kfun = estDeltaKdir}.
#'
#'A list of quadrats as required for argument \code{qsets} can be obtained by
#'function \code{\link{quadshilo}}.
#'}
#'\subsection{Details on the return value}{
#'The test returns an object belonging to classes \code{sostest} and \code{htest},
#'a list containing the following components:
#'\tabular{ll}{
#' \cr\code{statistic}\tab{the value of the test statistic,}
#' \cr\code{p.value}\tab{the p-value of the test,}
#' \cr\code{alternative}\tab{a character string describing the alternative hypothesis,}
#' \cr\code{method}\tab{a character string indicating what type of test was performed,}
#' \cr\code{data.name}\tab{a character string giving the name(s) of the data.}
#' \cr\code{Ksamples}\tab{a list of \code{\link{fdsample}}-objects, with elements 
#' \code{x}, \code{y} and \code{theo}}
#' }
#'}
#'
#'@export
#'@author Ute Hahn,  \email{ute@@imf.au.dk}
#' Hahn, U. (2012) A studentized permutation test for the comparison of spatial point  patterns.
#' \emph{Journal of the American Statistical Association} \strong{107} (498), 754--764.
#'
#' Hahn, U. and Jensen, E. B. V. (2013)
#' Inhomogeneous spatial point processes with hidden second-order stationarity.
#' \emph{CSGB preprint} 2013-7.
#' \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#'@seealso function \code{\link{tL2.permtest}} from package {fdnonpar} is used 
#'as test engine, \code{\link{quadshilo}} is used for setting up quadrat samples.

twosample.K.test <- function (x, y,
                      Kfun = estK,
                      r0, rlen = 256, 
                      ...,
                      use.tbar = FALSE,
                      nperm = 1000) 
{
  if(!is.ppsample(x) || !is.ppsample(y)) 
    stop("expect point pattern samples as arguments")
  AnisTest <- identical(Kfun, estDeltaKdir)
  dataname <- paste("point pattern samples", deparse(substitute(x)), 
                    "and", deparse(substitute(y)))
  rr <- seq(0, r0, length.out=rlen)
  Kfunx <- lapply(x, Kfun, r = rr, ...)
  Ktheo <- extract.fdsample(Kfunx[1], "theo")
  Kx <- extract.fdsample(Kfunx)
  Ky <- extract.fdsample(lapply(y, Kfun, r = rr, ...))
  if (use.tbar)  {
    Kxr <- fdsample(rr[rr>0], (Kx$fvals / rr)[rr>0, ])
    Kyr <- fdsample(rr[rr>0], (Ky$fvals / rr)[rr>0, ])
    testerg <-tL2.permtest(Kxr, Kyr, nperm = nperm, use.tbar = TRUE)
  } else {
    testerg <-tL2.permtest(Kx, Ky, nperm = nperm, use.tbar = FALSE)
  }
  method <- c(paste("Two sample studentized permutation test of identical K-functions,"),
              ifelse(AnisTest, "directional version, using Delta K_dir",
                       "isotropic version, using K_0"),
              paste("test statistic: ", if(use.tbar) "Tbar,"
                      else "T,", "upper integration bound:",r0),
                testerg$method[2])
  alternative <- c(paste("not the same K-function"),
                     if(AnisTest) ",\nbut different kinds of anisotropy")
  testerg$method <- method
  testerg$alternative <- alternative
  testerg$data.name <- dataname
  testerg$Ksamples <- list(x = Kx, y = Ky, theo = Ktheo)
  class(testerg) <- c("Ktest", "htest")
  return(testerg)
}
