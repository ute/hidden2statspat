#' Estimate the (template) K-function
#'
#' Estimates the \eqn{K}-function or template \eqn{K}-function of a point process.
#'
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
#' @param rmax optional numeric value: maximal r at which K-function is evaluated.
#' @param correction a character vector giving the edge correction type, may be
#'   any subset of \code{"border"},  \code{"isotropic"}, \code{"translate"}, \code{"none"}.
#' @param normpower an integer between 0 and 2. If \code{normpower} > 0, the
#'  intensity is normalized, see the Details.
#' @param ... optional arguments passed to \code{\link{as.sostyppp}}
#' @param max.ls.r optional, upper limit for argument \eqn{r} if \code{type="s"}.
#'
#' @details
#' By contrast to \code{\link{K.est}}, which returns an object of class \code{fv}, 
#' the function \code{estK} returns a \code{\link{funsample}}, which is a multivariate
#' function.
#'
#' \code{estK} takes both homogeneous and inhomogeneous point processes
#' that are supposed to be hidden second-order stationary (H&J, 2013). It then estimates
#' the so called \emph{template} \eqn{K}-function, the definition of which depends
#' on the type of second-order stationarity.
#' If \code{type} is not given, the last type of second-order stationarity assigned
#' to the point pattern \code{X} is used to determine how the template \eqn{K}-function 
#' is estimated.
#' If \code{X} has no type of second-order stationarity, it is assumed to be homogeneous.
#'
#' If \code{type} is given, but does not match the type of \code{X}, the function
#' \code{\link{as.sostyppp}} is called with arguments \ldots to ensure the correct
#' hidden second-order information.
#'
#' If  \code{normpower} > 0, the intensity is renormalized, so that \code{\link{estK}} 
#' yields similar results as \code{spatstat:\link[spatstat]{Kinhom}}. 
#' The intensity values \eqn{\lambda} are then multiplied by
#' \deqn{c^{normpower/2}}{c^(normpower/2),} where
#' \deqn{c = area(W)/sum_i(1/\lambda(x_i))}{c = area(W)/sum[i](1/lambda(x[i])).}
#'
#' The hidden \eqn{K}-function for \strong{reweighted} s.o. stationary point processes
#' delivers the same result as  \code{spatstat:\link[spatstat]{Kinhom}}, up to a subtle
#' difference for the border correction, where \code{estK} does not use fast optimized code.
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
#' as approximations. Here \code{estK} is similar to \code{spatstat:\link[spatstat]{Kscaled}}
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
#' # Parameters for plotting
#' Kstyle = simplist(col = list(theo = "black", trans = "red", iso = "blue", border = "green"),
#'                   lty = list(theo = "dashed"))
#' # For small point patterns, estK returns step functions by default
#' plot(estK(runifpoint(25), correction = c("trans", "iso")), Kstyle)
#' 
#' # compare homogeneous version of K.est and spatstat's Kest
#' # bronzefilter data are not marked as hidden second-order stationary
#' # set some colours and line styles for plotting
#' spatstatstyle = simplist(
#'     col = list(theo = "blue", trans = "red", iso = "black", border = "green"),
#'     lty = list(theo = "dotdash", trans = "dashed", iso = "solid", border = "dotted"))
#' plot(spatstat::Kest(bronzefilter))
#' plot(estK(bronzefilter), spatstatstyle)
#'
#' # Evaluate as reweighted, with default intensity estimate,
#' # compare with spatstat's Kinhom
#' plot(spatstat::Kinhom(bronzefilter, correction = "iso"))
#' plot(estK(reweighted(bronzefilter), correction = "iso"), Kstyle)
#' # There is a subtle difference because spatstat uses intensity renormalisation
#' # by default. We can do that, too:
#' plot(estK(bronzefilter, type = "w", normpower = 1, correction = "iso"), Kstyle)
#'
#' # Evaluate a rescaled version of the bronzefilter data:
#' plot(estK(rescaled(bronzefilter)))
#'
#' # The last given type is the one that counts
#' plot(estK(retransformed(rescaled(bronzefilter), backtrafo="gradx")))


estK <- function (X,
                  type,
                  rmax = NULL,
                  correction = c("border", "isotropic", "translate"),
                  normpower = 0,
                  ...,
                  max.ls.r = 3.0)
{
  # verifyclass(X, "ppp")
  npts <- npoints(X)
  stopifnot(npts > 1)

  if (missing(type))
  {
    if(is.sostyppp(X)) sostype <- currentType(X) else sostype <- NULL
    if (length(sostype) == 0)
      {
         X <- homogeneous(X)
         sostype <- "h"
      }
  }
  else {
    sostype <- type[1]
    X <- as.sostyppp(X, type = sostype, ...)
  }

  sostinfo <- attr(X, "sostinfo") 
  marx <- sostinfo$tmarks
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
    sostinfo <- attr(X, "sostinfo")
    marx <- sostinfo$tmarks
  }

  W <- X$window
  area <- area.owin(W)


  # get arguments r for K
  if (is.null(rmax))
  {
  if (scaling) {
    if(is.null(rmax)) rmax <- max.ls.r
    # calculate in the real world, without scaling
    rmaxdefault <- rmax.rule("K", W)
    # now transfer this to the scaled world
    rmaxdefault <- rmaxdefault * max(marx$invscale)
  }
  else {
    if(sostype == "w") lamax <- max(marx$intens) else lamax <- npts / area
    rmaxdefault <- rmax.rule("K", W, lamax)
    if(is.null(rmax)) rmax <- rmaxdefault
  }
  rmax <- rmaxdefault
  }  
  correction.given <- !missing(correction) && !is.null(correction)
  if (is.null(correction))
        correction <- c("border", "isotropic", "translate")
  correction <- match.arg(correction, c("border", "isotropic", "translate"),
                  several.ok = TRUE)
  correction <- implemented.for.K(correction, W$type, correction.given)


  # name of the function and its estimates
  typename <- sostype
  if (sostype == "hs") typename <- "s"
  # if (sostype == "h") typename <- "w"

  if (sostype %in% c("w", "s", "t"))
  {
    Kname <- paste("K[0]^{(", typename, ")}", sep="")
    Khatrname <- as.expression(substitute(widehat(K)[0]^(name)*(r), list(name = typename)))
    Ktheolab <- as.expression(substitute(K[0]^(name)*(r), list(name = typename)))
 # <<<< CHANGED HERE FOR PAPER >>>>
#    Ktheolab <- substitute(widehat(K)[0]^(name)*(r), list(name = typename))
  }
  else if (sostype == "h") {
    Kname <- "K"
    Khatrname <- expression(widehat(K)*(r))
    Ktheolab <- "K(r)"
  }
  else if (sostype == "hs") {
    Kname <- "K^*"
    Khatrname <- quote(widehat(K)^symbol('*')*(r))
    Ktheolab <- expression(K^symbol("*")*(r))
  }
 # Krname <- paste(Kname, "*(r)", sep="")
#  Khatrname <- paste(Khatname, "*(r)", sep="")
  legende <- function(AA = Khatrname, BB = "iso")
     as.expression(do.call("substitute", list(expression(AA*', '*scriptstyle(BB))[[1]], 
         list(AA=AA, BB = BB))))
  
# start funsample with CSR
  desc <- c("distance argument r", "theoretical Poisson %s")
  CSRlab <- if(typename == "s") "paste(K^symbol('*')*(r),', ',scriptstyle(CSR))"
            else  "paste(K*(r),', ',scriptstyle(CSR))"
  Ktheo <- urfunction(function(r) r^2 * pi, from = 0, to = rmax,
            xlab = "r", ylab = Ktheolab, legendtxt = CSRlab, main = "")
  KK <- list(theo = Ktheo)
  
  # identify all close pairs
 
  if(scaling) close <- lsclosepairs(X, rmax, invscale = marx$invscale)
  else close <- lsclosepairs(X, rmax)

  dIJ <- close$d
  eudIJ <- close$eud
  # compute weights for these pairs
  I <- close$i
  J <- close$j
  XI <- X[I]
  XJ <- X[J]
  
  rvals <- c(0, unique(sort(dIJ)))
  rsteps <- (length(rvals) > 1000) 
  if (rsteps) 
    rvals <- seq(0, rmax, length.out=501) 
    
# intensityweights. Intensities are equal to one if we deal with scaled processes.
# we use them here in the homogeneous / scaled case, too, and implement implicitely that
# infamous Poisson lambda^2 estimator n*(n-1)/area^2

  if (weighted) wIJ <- 1 / (marx$intens[J] * marx$intens[I])
  else wIJ <- 1 / marx$intens[J] * area / (npts - 1)

  Kmaker <- function(wh, rvals, cname = "") {
     Kst <- c(0,cumsum(wh) / area)
    if (rsteps) 
      Kfu <- approxfun (rvals, Kst) 
    else
      Kfu <- stepfun(rvals, c(0,Kst)) 
    Kfu <-  urfunction(Kfu, 
      xlim = c(0, rmax),  xlab = "r", ylab = Ktheolab, 
      legendtxt =  legende(BB = cname), do.points = FALSE,
      main = "")
    Kfu
  }

  if (any(correction == "none")) {
#  TODO: whist not really needed if rsteps is false
    wh <- whist(dIJ, rvals, wIJ)
    Kun <- Kmaker(wh, rvals, "none")
  }  else {
    Kun <- NULL
  }

  if (any(correction == "border")) {
    b <- bdist.points(X)
    if (scaling) b <- b * marx$invscale
    newwIJ <- wIJ
    if(scaling) newwIJ <- newwIJ * npts / area
    # horvitz-thompson weighting according to distance of point I to the boundary
    wtb <- npts / sapply(b, function(x) sum(b>=x))
    weightfun <- approxfun(b, wtb)
    HTweightsIJ <- weightfun(dIJ)*newwIJ
    useIJ <- dIJ <= b[I]
   
    wh <- whist(dIJ[useIJ], rvals, HTweightsIJ[useIJ])
    Kbord <- Kmaker(wh, rvals, "bord")
  } else {
    Kbord <- NULL
  }

  if (any(correction == "translate")) {
    edgewt <- edge.Trans(XI, XJ, paired = TRUE)
    totalwt <- edgewt * wIJ
    wh <- whist(dIJ, rvals, totalwt)
    Ktrans <-  Kmaker(wh, rvals, "trans")
  } else {
    Ktrans <- NULL
  }  

  if (any(correction == "isotropic")) {
    edgewt <- edge.Ripley(XI, matrix(eudIJ, ncol = 1))
    totalwt <- edgewt * wIJ
    wh <- whist(dIJ, rvals, totalwt)
    Kiso <- Kmaker(wh, rvals, "iso")
  } else {
    Kiso <- NULL
  } 
    
#  attr(K, "sostype") <- sostype
#  if (typename == "s") unitname (K) = NULL else unitname(K) <- unitname(X)
  allK <- list(theo = Ktheo, trans = Ktrans, iso = Kiso, border = Kbord)
  allK <- allK[!sapply(allK, is.null)]
  KK <- funsample(allK, 
    arglim = c(0, rmax),
    xlab = "r", ylab = Ktheolab, main = "" ) # TODO: sensible main
  attr(KK, "sostype") <- sostype
  return(KK)
}

#' Template L-function of a hidden 2nd-order stationary process
#'
#' Estimates the template \eqn{L}-function of a point process by transforming
#' the template \eqn{K}-function.
#'
#' @param ... Arguments for \code{\link{estK}}.
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
#' @seealso \code{\link{estK}}

estL <- function(...) {
  K <- estK(...)
  sostype <- attr(K, "sostype")
  # change here once there is a method for funsamples
  Kfuns <- attr(K, "funs")
  Lfuns <- vector("list", length(Kfuns))
  names (Lfuns) <- names(Kfuns)
  
  opt <- attr(K, "options")
  Klab <- opt$ylab
  if (is.character(Klab)) 
    Ltheolab <- sub("K", "L", Klab)
  else if (is.expression(Klab))
    Ltheolab <- as.expression(do.call("substitute", 
                          list(Klab[[1]], list(K="L"))))
  else Ltheolab <- Klab
  opt$ylab <- Ltheolab
  for (i in seq_along(Kfuns)) {
    optio <- attr(Kfuns[[i]], "options")
    optio$ylab <- Ltheolab
    if (is.stepfun(Kfuns[[i]])){
      rr <- knots(Kfuns[[i]])
      ll <- sqrt(Kfuns[[i]](rr)/pi)
      Lfuns[[i]] <- urfunction(stepfun(rr, c(0,ll)), optio)
    } else {
      if (names(Lfuns)[i] == "theo")
        Lfuns[[i]] <- urfunction(function(r) r, optio)
      else {
        # need to create a new object, since R might look up the wrong code.  
        # do this somehow q&d- no success with simple as.call and stuff things
        # TODO: make number of data points a package option
        argl <- attr(K, "arglim")
        rr <- seq(argl[1], argl[2], length.out = 501)
        ll <- sqrt(Kfuns[[i]](rr)/pi)
        Lfuns[[i]] <- urfunction(approxfun(rr, ll), optio)
      }
    }
  }
  
  LL <- funsample(Lfuns,
    arglim = attr(K, "arglim"),
    opt) 
  
  attr(LL, "sostype") <- sostype
  return(LL)
}
