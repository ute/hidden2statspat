# test on second-order stationarity

#'@title Test of second-order stationarity
#'@description Perform a permutation test of (hidden) second-order stationarity as 
#'described in Hahn & Jensen (2013). Returns an object of class \code{Kpermutest}
#'that also can be plotted, see the details.
#'@param x the point pattern to be tested, an object of class \code{sostyppp}
#'@param qsets optional subsamples to be compared, given as \code{list} 
#'    of (two) lists of  test quadrats or tessellations (objects of type 
#'    \code{\link{tess}}),  see Details
#'@param r0 numeric, the upper integration limit, see Details,
#'@param rlen optional, number of steps for numerical integration, defaults to 256,
#see Details,
#'@param Kfun optional \code{function}, the \eqn{K}-function to be used,
#'  either \code{\link{estK}} (default) or \code{\link{estDeltaKdir}},
#'@param ... optional parameters #for the function \code{\link{get.qsamples}}, and 
#'for the function \code{Kfun},# see Details
#'@param use.tbar logical, if true, a modified test statistic is used, see Details,
#'@param nperm number of random permutations, see Details.
#@param noTest optional logical, if \code{TRUE}, no test is run, only the point 
#pattern subsamples and \eqn{K}-function estimates are returned.
#'
#'@details The function tests if the point pattern \code{x} is (hidden) second-
#'order stationary. The pattern \code{x} is an object of type \code{sostyppp}, and
#'its type tag determines  which type of second-order stationarity will be tested.
#'\subsection{What the test does, and details on the parameters}{
#'The test is based on estimates of the K-function on two subsamples of the
#'pattern \code{x}. The two samples of estimated K-functions are compared by a 
#'permutation test. The test statistic is the integral over a squared Welch-t-statistic,
#'\deqn{T=\int_0^{r_0}\frac{(K_1(r)-K_2(r))^2}{s_1^2(r)/m_1 + 
#'   s_2^2(r)/m_2} d r}{T = integral [ (K_1(r)-K_2(r))^2 / (s_1^2(r)/m_1 + 
#'   s_2^2(r)/m_2)],} 
#'where \eqn{K_1(r)} and \eqn{K_2(r)} are the group means, and
#'\eqn{s_1^2(r), s_2^2(r)} are within group variances at a fixed argument \eqn{r}.
#'The integral spans an interval \eqn{[0, r_0]}. It is approximated by the mean
#'of the integrand over all \code{rlen} values of \eqn{r}, multiplied by the length
#'of the integral, i.e., by \code{r0}.
#'
#'A variant of the test statistic, \eqn{\bar T}{Tbar}, replaces the denominator in
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
#'To test against differences in local anisotropy, in particular to test the
#'hypotheses of locally rescaled and retransformed second-order stationarity 
#'against each other, H\&J (2013) propose to use the \eqn{\Delta K_{dir}}-function 
#'instead of the isotropic \eqn{K}-function. 
#'For this variant of the test, let \code{Kfun = estDeltaKdir}.
#'
#'A list of quadrats as required for argument \code{qsets} can be obtained by
#'function \code{\link{quadsets.hilo}}.
#'}
#\subsection{Specifiying test quadrats}{
#Samples of quadrats to be compared can be given explicitely, via the argument 
#\code{qsamples}. For best power, one should chose the two samples such that 
#one sample contains quadrats with high intensity, the other sample contains 
#quadrats with low intensity. Consequently, the list \code{qsamples} contains 
#two elements \code{hi} and \code{lo}. 
#These are
#either objects of ({spatstat}) class \code{\link{tess}}, or \code{list}s of
#objects of ({spatstat}-) class \code{\link{owin}}.
#
#If \code{qsamples} is not given, \code{sos.test} passes the optional \ldots 
#arguments to function \code{\link{get.qsamples}}. This function optionally 
#constructs quadrats, and divides them into sets with high and low intensity,
#conditioned on a minimum number \code{minpoints} of points. The arguments used
#are 
#\itemize{
#\item \code{quadrats} : pre-specified quadrats, a {spatstat}-\code{\link{tess}} object, 
#or a \code{list}s of objects of {spatstat}- \code{\link{owin}}s.
#\item\code{nx, ny, xbreaks, ybreaks, grad} : arguments for setting up quadrats,
#\item\code{minpoints} : the minimum number of points required (default: 20),
#}
#for details see \code{\link{get.qsamples}.
#}
#}  
#'\subsection{Details on the return value}{
#'The test returns an object belonging to classes \code{sostest} and \code{htest},
#'a list containing the following components:
#'\tabular{ll}{
#' \cr\code{statistic}\tab{the value of the test statistic,}
#' \cr\code{p.value}\tab{the p-value of the test,}
#' \cr\code{alternative}\tab{a character string describing the alternative hypothesis,}
#' \cr\code{method}\tab{a character string indicating what type of test was performed,}
#' \cr\code{data.name}\tab{a character string giving the name(s) of the data.}
#' \cr\code{ppsamples}\tab{a list of \code{{ppsample}}-objects, with elements 
#' \code{hi}, \code{lo} and \code{unused}}
#' \cr\code{Ksamples}\tab{a list of \code{\link{fdsample}}-objects, with elements 
#' \code{hi} and \code{lo}}
#' }
#'}
#'
#'@export
#'@author Ute Hahn,  \email{ute@@imf.au.dk}
#'@references
#'    Hahn, U. and Jensen, E. B. V. (2013)
#'    Inhomogeneous spatial point processes with hidden second-order stationarity.
#'    \emph{CSGB preprint} 2013-7.
#'    \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#'
#'@seealso function \code{\link{tL2.permtest}} from package {fdnonpar} is used 
#'as test engine, \code{\link{get.qsamples}} is used for setting up quadrat samples.
#'The quadrat subsamples or the $K$-function estimates can be plotted, see 
#'\code{\link{plot.sostest}}

sos.test <- function (x,
                      qsamples = NULL,
                      ...,
                      #quadrats = NULL,
                      #nx = NULL, ny = NULL, xbreaks = NULL, ybreaks = NULL,
                       #grad = "", minpoints = 20,
                      Kfun = estK,
                      r0 = NULL, rlen = 256, 
                      use.tbar = FALSE,
                      nperm = 25000,
                      noTest = FALSE) {
  NULL
}