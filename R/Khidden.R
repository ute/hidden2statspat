# the hidden K-function

#' Template L-function of a hidden 2nd-order stationary process
#'
#' Estimates the template \eqn{L}-function of a point process by transforming
#' the template \eqn{K}-function.
#'
#' @param \ldots Arguments for \code{\link{Khidden}}.
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
#' @seealso \code{\link{Khidden}}

Lhidden <- function(...) {
  K <- Khidden(...)
  htype <- attr(K, "htype")
  L <- eval.fv(sqrt(pmax.int(K, 0) / pi))
  # relabel the fv object
  L <- rebadge.fv(L, substitute(L[0]^(s)*(r), NULL), "widehat(L)[0]^(s)",
                  names(K), new.labl=attr(K, "labl"))
  #  L <- tweak.fv.entry(L, "theo", new.labl="paste(L[0]^{(s)}*(r),', ',scriptstyle(CSR))")
  L <- tweak.fv.entry(L, "theo", new.labl="paste(L[0]*(r),', ',scriptstyle(CSR))")
  attr(L, "htype") <-  htype
  return(L)
}

#' Locally rescaled K-function
#'
#' Estimates the locally rescaled \eqn{K}-function of a point process.
#' Modified (and frozen) version of \code{\link{Kest}}.
#'
#' @param X a point pattern, object of class \code{"ppp"}.
#' @param type character, the type of second-stationarity assumed: 
#'   \itemize{
#'        \item \code{"w"} reweighted
#'        \item \code{"t"} retransformed
#'        \item \code{"s"} locally rescaled
#'        \item \code{"h"} homogeneous, i.e. first order stationary
#'   }
#'   Only the first match is used.
#' @param r optional: vector of argument values \eqn{r} at which \eqn{K(r)} should
#'   be evaluated.
#' @param correction a character vector giving the edge correction type, may be
#'   any subset of \code{"border"},  \code{"isotropic"}, \code{"translate"}, \code{"none"}.
# @param r.max upper limit for argument \eqn{r}. The default 10 is huge - note
#  that the argument of the locally rescaled \eqn{K}-function corresponds to
#  \eqn{\r/\sqrt{\lambda}} in the non scaled case.
#'
#' @details
#'
#' For locally rescaled s.o. stationarity, all edge corrections are implemented 
#' as approximations. The translational edge
#' correction suffers from a small intrinsic bias in some cases of the locally
#' scaled point processes, depending on intensity and window shape. It is
#' recommended to use Ripley's isotropic edge correction, if possible.
#'
#' For more details, see \code{\link{Kscaled} and }
#' @export
#' @seealso \code{\link{Kest}}, \code{\link{Kinhom}}, \code{\link{Kscaled}}



Khidden <- function (X, 
                     type = c("w", "t", "s", "h", "hs"), 
                     r = NULL,
                     correction = c("border", "isotropic", "Ripley", "translate"),
                     max.ls.r = 6.0) 
{
  verifyclass(X, "ppp")
  
  matching <- matchtype(X, type)
  htype <- matching$htype
  marx <- matching$marx
  # algorithms to be used
  scaling <- htype %in% c("s", "hs")
  homogen <- htype %in% c("t", "h")
  weighted <- htype == "w"
  
  # name of the function and its estimates
  nametype = "htype"
  # if (htype == "hs") nametype ="s"
  # if (htype == "h") nametype ="w"
  Kname <- paste("K[0]^{(", nametype, ")}", sep="")
  Krname <- paste("K[0]^{(", nametype, ")}(r)", sep="")
  Khatname <- paste("widehat(K)[0]^(", nametype, ")", sep="")
  
    
  npts <- npoints(X)
  W <- X$window
  area <- area.owin(W)
  
  if (htype == "t") 
  {
    X <- extractRetransformed (X)
    # make uniform lambdas
    marx$lambda <- npts / area
  }
  
  
  if (scaling) { 
    # transition from scaled to real world
    maxrescale <- (marx$invscale)
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
    breaks <- handle.r.b.args(suggested.absr, NULL, W, rmaxdefault=rmaxdefault)
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
    if(htype == "w") lamax <- max(marx$lambda) else lamax <- npts / area
    rmaxdefault <- rmax.rule("K", W, lamax)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault = rmaxdefault)
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
                                                       translation = "translate", best = "best"), multi = TRUE)
  best.wanted <- ("best" %in% correction)
  correction <- implemented.for.K(correction, W$type, correction.given)
  
  demand.best <- correction.given && best.wanted

  
 
  
  
  Kdf <- data.frame(r = r, theo = pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  denom <- lambda2 * area
  K <- ratfv(Kdf, NULL, denom, "r", quote(K(r)), "theo", NULL, 
             alim, c("r", "%s[pois](r)"), desc, fname = "K", ratio = ratio)
  rmax <- max(r)
  close <- closepairs(X, rmax)
  DIJ <- close$d
  XI <- ppp(close$xi, close$yi, window = W, check = FALSE)
  if (any(correction == "none")) {
    wh <- whist(DIJ, breaks$val)
    numKun <- cumsum(wh)
    denKun <- lambda2 * area
    K <- bind.ratfv(K, data.frame(un = numKun), denKun, "hat(%s)[un](r)", 
                    "uncorrected estimate of %s", "un", ratio = ratio)
  }
  if (any(correction == "border" | correction == "bord.modif")) {
    b <- bdist.points(X)
    I <- close$i
    bI <- b[I]
    RS <- Kount(DIJ, bI, b, breaks)
    if (any(correction == "bord.modif")) {
      denom.area <- eroded.areas(W, r)
      numKbm <- RS$numerator
      denKbm <- lambda2 * denom.area
      K <- bind.ratfv(K, data.frame(bord.modif = numKbm), 
                      data.frame(bord.modif = denKbm), "hat(%s)[bordm](r)", 
                      "modified border-corrected estimate of %s", "bord.modif", 
                      ratio = ratio)
    }
    if (any(correction == "border")) {
      numKb <- RS$numerator
      denKb <- lambda * RS$denom.count
      K <- bind.ratfv(K, data.frame(border = numKb), data.frame(border = denKb), 
                      "hat(%s)[bord](r)", "border-corrected estimate of %s", 
                      "border", ratio = ratio)
    }
  }
  if (any(correction == "translate")) {
    XJ <- ppp(close$xj, close$yj, window = W, check = FALSE)
    edgewt <- edge.Trans(XI, XJ, paired = TRUE)
    wh <- whist(DIJ, breaks$val, edgewt)
    numKtrans <- cumsum(wh)
    denKtrans <- lambda2 * area
    h <- diameter(as.rectangle(W))/2
    numKtrans[r >= h] <- NA
    K <- bind.ratfv(K, data.frame(trans = numKtrans), denKtrans, 
                    "hat(%s)[trans](r)", "translation-corrected estimate of %s", 
                    "trans", ratio = ratio)
  }
  if (any(correction == "isotropic")) {
    edgewt <- edge.Ripley(XI, matrix(DIJ, ncol = 1))
    wh <- whist(DIJ, breaks$val, edgewt)
    numKiso <- cumsum(wh)
    denKiso <- lambda2 * area
    h <- diameter(W)/2
    numKiso[r >= h] <- NA
    K <- bind.ratfv(K, data.frame(iso = numKiso), denKiso, 
                    "hat(%s)[iso](r)", "Ripley isotropic correction estimate of %s", 
                    "iso", ratio = ratio)
  }
  
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  nama <- nama[!(nama %in% c("r", "rip", "ls"))]
  fvnames(K, ".") <- nama
  unitname(K) <- unitname(X)
  
  return(K)
}



# get r-argument ----------------------------------------------------------

# replacement for handle.r.b.args and rmax.rule for locally scaled point processes
# @param r 
# @param type assumed type of second-order stationarity
# @return a list (\code{htype}, \code{marx}) with matched type and mark data frame of \code{X}
#' @rdname sosspp-internal
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
