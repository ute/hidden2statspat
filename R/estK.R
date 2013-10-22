# the hidden K-function

#' Template L-function of a hidden 2nd-order stationary process
#'
#' Estimates the template \eqn{L}-function of a point process by transforming
#' the template \eqn{K}-function.
#'
#' @param ... Arguments for \code{\link{K.est}}.
#' @details
#' For a reweighted or retransformed second-order stationary Poisson point process,
#' the theoretical value of the template \eqn{L}-function is \deqn{L_0(r)=r}{L_0(r)=r}.
#'
#' Under locally rescaled second-order stationarity,
#' this holds only approximately for large arguments, if the intensity varies strongly.
#' Furthermore, approximations when calculating
#' locally scaled distances and simplifications in the estimation of
#' the locally scaled \eqn{K}-function may add a bias to the estimates of the locally scaled
#' \eqn{L}-function for inhomogeneous Poisson point processes.
#' The bias depends upon how fast the intensity function varies. Usually, the approximation
#' holds very well for reasonably small arguments.
#' @export
#' @seealso \code{\link{K.est}}

L.est <- function(...) {
  K <- K.est(...)
  sostype <- attr(K, "sostype")
  L <- eval.fv(sqrt(pmax.int(K, 0) / pi))
  attr(L, "sostype") <-  sostype
  if (sostype %in% c("w", "s", "t"))
  {
    Lhatname <- paste("widehat(L)[0]^(", sostype, ")", sep="")
    Ltheolab <- substitute(L[0]^(name)*(r), list(name = sostype))
  }
  else if (sostype == "h") {
    Lhatname <- "widehat(L)"
    Ltheolab <- "L(r)"
  }
  else if (sostype == "hs") {
    Lhatname <- "widehat(L)^symbol('*')"
    Ltheolab <- expression(L^symbol('*')*(r))
  }

  # relabel the fv object
  L <- rebadge.fv(L, Ltheolab, Lhatname, names(K), new.labl=attr(K, "labl"))
  CSRlab <- if(sostype %in% c("s", "hs")) "paste(L^symbol('*')*(r),', ',scriptstyle(CSR))"
             else  "paste(L*(r),', ',scriptstyle(CSR))"
  L <- tweak.fv.entry(L, "theo", new.labl=CSRlab)

  return(L)
}

#' Estimate the (template) K-function
#'
#' Estimates the \eqn{K}-function or template \eqn{K}-function of a point process.
#' Modified (and frozen) version of \code{spatstat:\link[spatstat]{Kest}}.
#'
#' @param X a point pattern, object of class \code{"ppp"}.
#' @param type optional character, the type of second-order stationarity assumed:
#'   \itemize{
#'        \item \code{"w"} reweighted
#'        \item \code{"t"} retransformed
#'        \item \code{"s"} locally rescaled
#'        \item \code{"h"} homogeneous, i.e. first order stationary
#'        \item \code{"hs"} homogeneous, but evaluated as scaled \eqn{K}-function
#'   }
#'   Only the first match is used.
#' @param r optional: vector of argument values \eqn{r} at which \eqn{K(r)} should
#'   be evaluated.
#' @param correction a character vector giving the edge correction type, may be
#'   any subset of \code{"border"},  \code{"isotropic"}, \code{"translate"}, \code{"none"}.
#' @param normpower an integer between 0 and 2. If \code{normpower} > 0, the
#'  intensity is normalized, see the Details.
#' @param ... optional arguments passed to \code{\link{as.sostyppp}}
#' @param max.ls.r optional, upper limit for argument \eqn{r} if \code{type="s"}.
#'
#' @details
#' If applied to homogeneous point patterns, \code{sostatpp}'s function \code{K.est}
#' is a simplified version of \code{spatstat:\link[spatstat]{Kest}}.
#' In particular, it does not provide any variance estimates, and does not choose
#' the edge correction according to the number of points in the data set. See the excellent
#' man page of \code{spatstat:\link[spatstat]{Kest}} for information on the \eqn{K}-function
#' and edge correction.
#'
#' The present version of \code{K.est} also applies to inhomogeneous point processes
#' that are supposed to be hidden second-order stationary (H&J, 2013). It then estimates
#' the so called \emph{template} \eqn{K}-function, the definition of which depends
#' on the type of second-order stationarity.
#' If \code{type} is not given, the last type of second-order stationarity assigned
#' to \code{X} is used to determine how the template \eqn{K}-function is estimated.
#' If \code{X} has no type of second-order stationarity, it is assumed to be homogeneous.
#'
#' If \code{type} is given, but does not match the type of \code{X}, the function
#' \code{\link{as.sostyppp}} is called with arguments ... to ensure the correct
#' hidden second-order information.
#'
#' If  \code{normpower} > 0, the intensity is renormalized, so that \code{\link{K.est}} yields similar results as
#' \code{spatstat:\link[spatstat]{Kinhom}}. The intensity values \eqn{\lambda} are then multiplied by
#' \deqn{c^{normpower/2}}{c^(normpower/2),} where
#' \deqn{c = area(W)/sum_i(1/\lambda(x_i))}{c = area(W)/sum[i](1/lambda(x[i])).}
#'
#' The hidden \eqn{K}-function for \strong{reweighted} s.o. stationary point processes
#' delivers the same result as  \code{spatstat:\link[spatstat]{Kinhom}}, up to a subtle
#' difference for the border correction: \code{K.est} does not use fast optimized code.
#'
#' If \code{X} is typed \strong{retransformed} s.o. stationary, the \eqn{K}-function of
#' the \code{\link{backtransformed}} point pattern is returned, which corresponds to
#' {spatstat}'s function \code{spatstat:\link[spatstat]{Kest}}
#'
#' For \strong{locally rescaled} s.o. stationarity, the locally scaled interpoint
#' distances are computed using an approximation proposed by Hahn (2007):
#' The Euclidean distance between two points is multiplied by the
#' average of the square roots of the intensity values at the two points.
#' Similarly, all edge corrections are implemented
#' as approximations. Here \code{K.est} is similar to \code{spatstat:\link[spatstat]{Kscaled}}
#'  The translational edge
#' correction suffers from a small intrinsic bias in some cases of locally
#' rescaled s.o.s., depending on intensity and window shape. It is
#' recommended to use Ripley's isotropic edge correction, if possible.
#'
#' Note that the argument of the locally rescaled \eqn{K}-function corresponds to
#'  \deqn{r/\sqrt{\lambda}}{r / sqrt(lambda)} in the non scaled case.
#'  The default maximum \code{max.ls.r} = 3 is thus quite large.
#'
#' @export
#' @seealso the \pkg{spatstat} functions \code{\link[spatstat]{Kest}}, \code{\link[spatstat]{Kinhom}}, \code{\link[spatstat]{Kscaled}}
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @references
#'    Prokesova, M., Hahn, U. and Vedel Jensen, E.B. (2006)
#'    Statistics for locally scaled point patterns.
#'    In A. Baddeley, P. Gregori, J. Mateu, R. Stoica and D. Stoyan (eds.)
#'    \emph{Case Studies in Spatial Point Pattern Modelling.}
#'    Lecture Notes in Statistics 185. New York: Springer Verlag. Pages 99--123
#'
#'    Hahn, U. (2007) \emph{Global and Local Scaling in the Statistics of
#'    Spatial Point Processes}. Habilitationsschrift, Universitaet Augsburg.
#'
#'    Hahn, U. and Jensen, E. B. V. (2013)
#'    Inhomogeneous spatial point processes with hidden second-order stationarity.
#'    \emph{CSGB preprint} 2013-7.
#'    \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#'
#' @examples
#' # compare homogeneous version of K.est and spatstat's Kest
#' # bronzefilter data are not marked as hidden second-order stationary
#' plot(spatstat::Kest(bronzefilter))
#' plot(K.est(bronzefilter))
#'
#' # Evaluate as reweighted, with default intensity estimate,
#' # compare with spatstat's Kinhom
#' plot(spatstat::Kinhom(bronzefilter))
#' plot(K.est(reweighted(bronzefilter)))
#' # There is a subtle difference because spatstat uses intensity renormalisation
#' # by default. We can do that, too:
#' plot(K.est(bronzefilter, type = "w", normpower = 1))
#'
#' # Evaluate a rescaled version of the bronzefilter data:
#' plot(K.est(rescaled(bronzefilter)))
#'
#' # The last given type is the one that counts
#' plot(K.est(retransformed(rescaled(bronzefilter), backtrafo="gradx")))


K.est <- function (X,
                  type,
                  r = NULL,
                  correction = c("border", "isotropic", "Ripley", "translate"),
                  normpower = 0,
                  ...,
                  max.ls.r = 3.0)
{
  # verifyclass(X, "ppp")
  npts <- npoints(X)
  stopifnot(npts > 1)

  if (missing(type))
  {
    if(is.sostyppp(X)) sostype <- currenttype(X) else sostype <- NULL
    if (length(sostype) == 0)
      {
         X <- ashomogeneous(X)
         sostype <- "h"
      }
  }
  else {
    sostype <- type[1]
    X <- as.sostyppp(X, type = sostype, ...)
  }

  marx <- X$sostinfo$tmarks
  if (normpower != 0) {
    stopifnot ((1 <= normpower) & (normpower <= 2))
    if (!is.null(marx$intens)){
      renorm.factor <-  (sum(1 / marx$intens) / (area.owin(X)))^(normpower / 2)
      marx$intens <- marx$intens * renorm.factor
    }
    if(!is.null(marx$invscale)){
      renorm.factor <-  (sum(1 / marx$invscale^2) / (area.owin(X)))^(normpower / 4)
      marx$invscale <- marx$invscale * renorm.factor
    }
  }

  # algorithms to be used
  scaling <- sostype %in% c("s", "hs")
  homogen <- sostype %in% c("t", "h")
  weighted <- sostype == "w"

  # modify point pattern to "become" the template, if retransformed or rescaled

  if (scaling) marx$intens <- rep(1, npts) # refers to unit rate template

  if (sostype == "t")
  {
    X <- backtransformed (X)
    # make uniform intensities
    # marx$intens <- rep(npts / area, npts)
    marx <- X$sostinfo$tmarks
  }

  W <- X$window
  area <- area.owin(W)


  # get arguments r for K

  if (scaling) {
    # transition from scaled to real world
    maxrescale <- max(marx$invscale)
    if(!is.null(r)) {
      suggested.absr <- r/maxrescale
      r.max <- max(r)
    }
    else {
      suggested.absr <-NULL
      r.max <- max.ls.r
    }
    # calculate in the real world, without scaling
    rmaxdefault <- rmax.rule("K", W)
    breaks <- NULL
    breaks <- handle.r.b.args(suggested.absr, breaks, W, rmaxdefault=rmaxdefault)
    # get back to scaled world
    rmaxdefault <- rmaxdefault * maxrescale
    breaks$val <- breaks$val * maxrescale
    breaks$r <- breaks$r * maxrescale
    r <- breaks$r
    rmax <- max(r)
    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault, max.ls.r))
  }
  else {
    if(sostype == "w") lamax <- max(marx$intens) else lamax <- npts / area
    rmaxdefault <- rmax.rule("K", W, lamax)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max
    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
  }


  correction.given <- !missing(correction) && !is.null(correction)
  if (is.null(correction))
        correction <- c("border", "isotropic", "Ripley", "translate")
  correction <- pickoption("correction", correction, c(none = "none",
                                                       border = "border", bord.modif = "bord.modif", isotropic = "isotropic",
                                                       Ripley = "isotropic", trans = "translate", translate = "translate",
                                                       translation = "translate"), multi = TRUE)
  correction <- implemented.for.K(correction, W$type, correction.given)


  # name of the function and its estimates
  typename = sostype
  if (sostype == "hs") typename ="s"
  # if (sostype == "h") typename ="w"

  if (sostype %in% c("w", "s", "t"))
  {
    Kname <- paste("K[0]^{(", typename, ")}", sep="")
    Khatname <- paste("widehat(K)[0]^(", typename, ")", sep="")
    Ktheolab <- substitute(K[0]^(name)*(r), list(name = typename))
  }
  else if (sostype == "h") {
    Kname <- "K"
    Khatname <- "widehat(K)"
    Ktheolab <- "K(r)"
  }
  else if (sostype == "hs") {
    Kname <- "K^*"
    Khatname <- "widehat(K)^symbol('*')"
    Ktheolab <- expression(K^symbol("*")*(r))
  }
  Krname <- paste(Kname, "*(r)", sep="")

  # start data frame with CSR
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", Ktheolab,
            "theo", , alim, c("r",Krname), desc, fname=Kname)
  CSRlab <- if(typename == "s") "paste(K^symbol('*')*(r),', ',scriptstyle(CSR))"
            else  "paste(K*(r),', ',scriptstyle(CSR))"
  K <- tweak.fv.entry(K, "theo", new.labl=CSRlab)


 # K <- tweak.fv.entry(K, "theo", new.labl="paste(Krname,', ',scriptstyle(CSR))")


  # identify all close pairs
  rmax <- max(r)

  if(scaling) close <- lsclosepairs(X, rmax, invscale = marx$invscale)
  else close <- lsclosepairs(X, rmax)

  dIJ <- close$d
  eudIJ <- close$eud
  # compute weights for these pairs
  I <- close$i
  J <- close$j
  XI <- X[I]
  XJ <- X[J]

# intensityweights. No worries, mate, intensities are one if we deal with scaled processes.
# we use them here in the homogeneous / scaled case, too, and implement implicitely that
# infamous Poisson lambda^2 estimator n*(n-1)/area^2

  if (weighted) wIJ <- 1 / (marx$intens[J] * marx$intens[I])
  else wIJ <- 1 / marx$intens[J] * area / (npts - 1)


  if (any(correction == "none")) {
    wh <- whist(dIJ, breaks$val, wIJ)
    Kun <- cumsum(wh) / area
    K <- bind.fv(K, data.frame(un = Kun), "%s[un](r)",
                    "uncorrected estimate of %s", "un")
  }

  if (any(correction == "border" | correction == "bord.modif")) {
    b <- bdist.points(X)
    if (scaling) b <- b * marx$invscale
    bI <- b[I]
    newwIJ <- wIJ
    if(scaling) newwIJ <- newwIJ * npts / area
    RS <- Kwtsum(dIJ, bI, newwIJ, b, w = 1/marx$intens, breaks)
    if (any(correction == "border")) {
      Kb <- RS$ratio
      K <- bind.fv(K, data.frame(border=Kb),"%s*(r)", # "%s[bord](r)",
                   "border-corrected estimate of %s",
                   "border")
      K <- rebadge.fv(K, Ktheolab, Khatname)
      K <- tweak.fv.entry(K, "border", new.labl="paste(%s*(r),', ',scriptstyle(bord))")
    }
    if (any(correction == "bord.modif")) {
      Kbm <- RS$numerator / eroded.areas(W, r)
      K <- bind.fv(K, data.frame(bord.modif = Kbm),"%s*(r)", # "%s[bordm](r)",
                   "modified border-corrected estimate of %s",
                   "bord.modif")
      K <- rebadge.fv(K, Ktheolab, Khatname)
      K <- tweak.fv.entry(K, "bord.modif", new.labl="paste(%s*(r),', ',scriptstyle(bord.mod))")
    }
  }

  if (any(correction == "translate")) {
    edgewt <- edge.Trans(XI, XJ, paired = TRUE)
    totalwt <- edgewt * wIJ
    wh <- whist(dIJ, breaks$val, totalwt)
    Ktrans <- cumsum(wh) / area
    K <- bind.fv(K, data.frame(trans=Ktrans), "%s[trans](r)",
                 "translation-corrected estimate of %s",
                 "trans")
    K <- rebadge.fv(K, Ktheolab, Khatname)
    K <- tweak.fv.entry(K, "trans", new.labl="paste(%s*(r),', ',scriptstyle(trans))")
  }

  if (any(correction == "isotropic")) {
    edgewt <- edge.Ripley(XI, matrix(eudIJ, ncol = 1))
    totalwt <- edgewt * wIJ
    wh <- whist(dIJ, breaks$val, totalwt)
    Kiso <- cumsum(wh) / area
    K <- bind.fv(K, data.frame(iso=Kiso), "%s[iso](r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso")
    K <- rebadge.fv(K, Ktheolab, Khatname)
    K <- tweak.fv.entry(K, "iso", new.labl="paste(%s*(r),', ',scriptstyle(iso))")
  }

  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  nama <- nama[nama != "r"] #!(nama %in% c("r", "rip", "ls"))]
  fvnames(K, ".") <- nama
  attr(K, "sostype") <- sostype
  if (typename == "s") unitname (K) = NULL else unitname(K) <- unitname(X)
  return(K)
}



# get r-argument ----------------------------------------------------------

# replacement for handle.r.b.args and rmax.rule for locally scaled point processes
# @param r
#' @rdname sostatpp-internal
#' @keywords internal
#'

ls.r.args <- function (r = NULL, eps = NULL, rmaxdefault = 4, minscale = NULL, maxscale)
{
  if (is.null(minscale)) minscale <- max
  if (!is.null(r) && !is.null(breaks))
    stop(paste("Do not specify both", sQuote("r"), "and",
               sQuote("breaks")))
  if (!is.null(breaks)) {
    breaks <- as.breakpts(breaks)
  }
  else if (!is.null(r)) {
    breaks <- breakpts.from.r(r)
  }
  else {
    rmax <- if (missing(rmaxdefault))
      diameter(as.rectangle(window))
    else rmaxdefault
    if (is.null(eps)) {
      if (!is.null(window$xstep))
        eps <- window$xstep/4
      else eps <- rmax/512
    }
    breaks <- make.even.breaks(rmax, bstep = eps)
  }

  # r <- breaks$r
  # rmax <- breaks$max

  return(breaks)
}
