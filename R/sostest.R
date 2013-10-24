# test of second-order stationarity

#'@title Test of second-order stationarity
#'@description Perform a permutation test of (hidden) second-order stationarity as
#'described in Hahn & Jensen (2013). Returns an object of class \code{Kpermutest}
#'that also can be plotted, see the details.
#'@param x the point pattern to be tested, an object of class \code{sostyppp}
#'@param qsets quadrat sets on which the point pattern is to be compared, given as
#'     \code{list}
#'    of (two) lists of  test quadrats or tessellations (objects of type
#'    \code{\link{tess}}). The list contains entries with names \code{hi}
#'    and \code{lo}, optionally also \code{unused}, as produced by  \code{\link{twoquadsets}}.
#'@param rmax numeric, the upper integration limit, see Details,
#'@param rlen optional, number of steps for numerical integration, defaults to 256,
#see Details,
#'@param Kfun optional \code{function}, the \eqn{K}-function to be used,
#'  either \code{\link{K.est}} (default) or \code{\link{DeltaKdir.est}},
#'@param ... optional parameters for function \code{Kfun}. To speed up
#'calculations, it is recommended to give an explicite \code{correction} argument.
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
#'\deqn{T=\int_0^{r_{max}}\frac{(K_1(r)-K_2(r))^2}{s_1^2(r)/m_1 +
#'   s_2^2(r)/m_2} d r}{T = integral [ (K_1(r)-K_2(r))^2 / (s_1^2(r)/m_1 +
#'   s_2^2(r)/m_2)],}
#'where \eqn{K_1(r)} and \eqn{K_2(r)} are the group means, and
#'\eqn{s_1^2(r), s_2^2(r)} are within group variances at a fixed argument \eqn{r}.
#'The integral spans an interval \eqn{[0, r_{max}]}. It is approximated by the mean
#'of the integrand over all \code{rlen} values of \eqn{r}, multiplied by the length
#'of the integral, i.e., by \code{rmax}.
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
#'To test against differences in local anisotropy, in particular to test the
#'hypotheses of locally rescaled and retransformed second-order stationarity
#'against each other, H\&J (2013) propose to use the \eqn{\Delta K_{dir}}-function
#'instead of the isotropic \eqn{K}-function.
#'For this variant of the test, let \code{Kfun = DeltaKdir.est}.
#'
#'A list of quadrats as required for argument \code{qsets} can be obtained by
#'function \code{\link{twoquadsets}}.
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
# \cr\code{ppsamples}\tab{a list of \code{{ppsample}}-objects, with elements
# \code{hi}, \code{lo} (and \code{unused}, if present)}
# takes muc too much space for some examples with intensity as im
#' \cr\code{Ksamples}\tab{a list of \code{\link{fdsample}}-objects, with elements
#' \code{hi}, \code{lo} and \code{theo}}
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
#'as test engine, \code{\link{twoquadsets}} is used for setting up quadrat samples.
#'The quadrat subsamples or the $K$-function estimates can be plotted, see
#'\code{\link{plot.Ktest}}

sos.test <- function (x,
                      qsets = NULL,
                      Kfun = K.est,
                      rmax, rlen = 256,
                      ...,
                      use.tbar = FALSE,
                      nperm = 1000) {
  if(!is.ppp(x)) stop("not a point pattern")
  AnisTest <- identical(Kfun, DeltaKdir.est)
  dataname <- paste( "point pattern",deparse(substitute(sost)))
  type <- currenttype(x)
  if (length(type) == 0) stop ("unknown type of hidden 2nd-order stationarity")
  typename <- .TYPENAMES[which(.TYPES == type)]

  stopifnot (is.list(qsets) && length(qsets) > 1)
  stopifnot (all(c("hi", "lo") %in% names(qsets)))
  pp.hi <- ppsubsample(x, qsets$hi)
  pp.lo <- ppsubsample(x, qsets$lo)
  pp.unused <- if (!is.null(qsets$unused)) ppsubsample(x, qsets$unused) else NULL

  testerg <- twosample.K.test(pp.hi, pp.lo, Kfun = Kfun, rmax = rmax, ...,
                              use.tbar = use.tbar, nperm = nperm)
  names(testerg$Ksamples) <- c("theo", "hi", "lo")
 
# 
#   testerg$ppsamples <- list(hi = pp.hi, lo = pp.lo, unused = pp.unused)

  testerg$data.name <- dataname

  testerg$method <- c(paste("Studentized permutation test of",
          typename ," hidden second-order stationarity,"),
          ifelse(AnisTest, "directional version, using Delta K_dir",
                "isotropic version, using K_0"),
          paste("test statistic: ", if(use.tbar) "Tbar,"
               else "T,", "upper integration bound: ", rmax),
           testerg$method[2] )
  testerg$alternative <- c(paste("not the same", typename, "K-function"),
                      if(AnisTest) ",\nbut different kinds of anisotropy")
  testerg
}


#'@title Plot K-functions used in test on second-order stationarity
#'@description Plot the K-function estimates used
#'in a test on second-order stationarity, \code{\link{sos.test}}.
#'@param x result of a \code{\link{sos.test}} or a \code{\link{twosample.K.test}},
#'an object of class \code{"Ktest"}
#'@param styles named list of \code{\link{style}} lists, defines how  the K-function
#'estimates are plotted, see Details,
#'@param theostyle plot style for the reference \eqn{K}-function of a Poisson point process.
#'To supress plotting of the reference curve, let \code{theostyle = NULL}.
#'@param ... further arguments passed to plot methods,
#'@param labline numeric, controls the position of the axis labels: first 
#'entry in graphic parameter \code{mgp}. 
#@param mean.thicker optional numeric, multiplier for the line width of the group mean functions,
#@param mean.alpha optional numeric, alpha value for the colour of the group mean functions.
#'@details
#If \code{plotquads} is \code{FALSE}, t
#'The estimates of the \eqn{K}-function on the quadrats are plotted together
#'with the group means. The \eqn{K}-function
#'of a Poisson point process is also present in the plot, unless \code{theostyle = NULL}.
#'
#If \code{plotquads} is \code{TRUE}, a plot of the point pattern, divided into
#the quadrats that were used in the test, is generated, with background colours
#according to \code{hilostyles}.
#
#'Colours and line attributes for plotting K-function estimates
#'are controlled by argument \code{styles}, which contains two elements.
#'For plotting the result of a \code{sos.test}, the have to be named \code{hi}
#'and \code{lo}, for the result of a \code{twosample.K.test}, the names
#'are \code{x} and \code{y}. The elements themselve are \code{\link{style}} lists
#'with elements
#'\tabular{ll}{
#'\code{col} \tab colour for the individual \eqn{K}-functions
#\cr\code{col.mea}\tab optional, colour for the mean \eqn{K}-function, defaults to \code{col}
#\cr\code{col.env}\tab optional, colour for the envelope, defaults to \code{col}
#\cr\code{col.win}\tab optional, colour for the quadrats, defaults to \code{col}
#\cr\code{col.pts}\tab optional, colour for the points in quadratsplot, defaults to
#\code{par("col")}
#'\cr\code{alpha} \tab optional, alpha-value for the individual \eqn{K}-functions,
#'defaults to \code{0.5}
#\cr\code{alpha.win} \tab optional, alpha-value for the quadrats, defaults to \code{0.5 * alpha}
#'\cr\code{lwd}\tab optional, line width for plot of \eqn{K}-functions
#\cr\code{lwd.sum}\tab optional line width for the mean \eqn{K}-function, defaults to \code{2 * lwd}
#'\cr\code{lty}\tab optional, line type for plot of \eqn{K}-functions
#\cr\code{lty.sum}\tab optional line type for the mean \eqn{K}-function, defaults to \code{lty}
#'}
#'
#'@export plot.Ktest
#'@method plot Ktest
#'@author Ute Hahn,  \email{ute@@imf.au.dk}
#'@examples
#'# testing beilschmiedia pattern on reweighted second-order stationarity
#'bei.ml <- reweighted(bei, intensity = bei.intens.maxlik)
#'bei.quads <- twoquadsets(bei, nx = 8, ny = 4, minpoints = 30)
#'beitest <- sos.test(bei.ml, qsets = bei.quads, rmax =25 )
#'beistyle <- list(hi = style(col = "red", alpha = .5), lo = style(col = "blue", alpha = .5))
#'
#'plot(beitest, beistyle, main = "bei.ml: K estimated on quadrats", ylim = c(0,3000))

plot.Ktest <- function(x, styles,
                      theostyle = style(lty = "dotted", col = "black", alpha = 1),
                       ..., labline = 2.4)
  #,
  #                     mean.thicker = 2, mean.alpha = 1)
{
  Ksamp <- x$Ksamples
  if (is.null(theostyle)) Ksamp$theo <- NULL
  if (missing(styles)) styles <- list()
  styles$theo <- theostyle
  # don't throw warnings because of extra arguments in the styles
  okstyles <- lapply(styles, matching, plot.fdsample, .graphparams)
  textlines <- par("mgp")
  textlines[1] <- labline
  par(mgp = textlines)

  allrange <- sapply(Ksamp, rangexy, finite = TRUE)

  dotargs <- style(...)
  if (is.null(dotargs$xlim))
    dotargs$xlim <- range(range(allrange["x",]))
  if (is.null(dotargs$ylim))
    dotargs$ylim <- range(range(allrange["y",]))

  if (!is.null(Ksamp$theo)) {
    if (is.null(theostyle$lwd)) {
      if(!is.null(dotargs$lwd)) theostyle$lwd <- 2*dotargs$lwd
      else theostyle$lwd <- 2*par("lwd")
    }
    splot(Ksamp$theo, dotargs, theostyle, add = FALSE)
    add = TRUE
  } else {
    add = FALSE
  }

  Ksamp$theo <- NULL
  lplot(Ksamp, styles, dotargs, add = add, allinone = TRUE, .plotmethod = "summaryplot")

}
