#
# Klocscaled.R Estimation of K function for locally-scaled process
#

#' Locally rescaled L-function
#' 
#' Estimates the locally rescaled template \eqn{L}-function of a point process by transforming
#' the locally scaled \eqn{K}-function.
#' 
#' @param \ldots Arguments for \code{\link{Klocscaled}}.
#' @details 
#' For a Poisson point process, the theoretical value of the locally rescaled 
#' \eqn{L}-function is \eqn{L_0(r)=r}. Due to approximations when calculating 
#' locally scaled distances and due to a simplifications in the estimation of 
#' the locally scaled \eqn{K}-function, the estimates of the locally scaled 
#' \eqn{L}-function may be biased for inhomogeneous Poisson point processes. 
#' The bias depends upon how fast the intensity function varies.
#' @export
#' @seealso \code{\link{Klocscaled}}

Llocscaled <- function(...) {
  K <- Klocscaled(...)
  L <- eval.fv(sqrt(pmax.int(K,0)/pi))
  # relabel the fv object
  L <- rebadge.fv(L, substitute(L[0]^(s)*(r), NULL), "widehat(L)[0]^(s)",
                  names(K), new.labl=attr(K, "labl"))
#  L <- tweak.fv.entry(L, "theo", new.labl="paste(L[0]^{(s)}*(r),', ',scriptstyle(CSR))") 
    L <- tweak.fv.entry(L, "theo", new.labl="paste(L[0]*(r),', ',scriptstyle(CSR))") 
return(L)  
}

#' Locally rescaled K-function
#' 
#' Estimates the locally rescaled \eqn{K}-function of a point process.
#' Modified, fixed (and frozen) version of \code{\link{Kscaled}}.
#' 
#' @param X a point pattern, object of class \code{"ppp"}.
#' @param lambda optional: intensity function: a vector of values at 
#'   the points of \code{X}, a pixel image, or a \code{function (x,y)}. If not 
#'   given, it is estimated by \code{\link{density.ppp}}.
#' @param \ldots arguments for \code{\link{density.ppp}} if \code{lambda} is omitted.
#' @param r optional: vector of argument values \eqn{r} at which \eqn{K(r)} should
#'   be evaluated.
#'  @param correction a character vector giving the edge correction type, may be 
#'   any of \code{"border"},  \code{"isotropic"}, \code{"translate"}, code{"none"}.
#' @param renormalise logical, defaults to FALSE. If \code{TRUE}, the intensity is renormalised, see `Details'.
#' @param normpower  integer (either 1 or 2). Normalisation power, see `Details'.
#' @param sigma,varcov optional arguments for \code{\link{density.ppp}}
#'  to control kernel bandwidth when \code{lambda} is estimated.
#' @param r.max upper limit for argument \eqn{r}. The default 10 is huge - note
#'  that the argument of the locally rescaled \eqn{K}-function corresponds to 
#'  \eqn{\r/\sqrt{\lambda}} in the non scaled case.
#'   
#' @details 
#' This version of the locally rescaled \eqn{K}-function has been extended by a 
#' renormalisation of the intensity, as in \code{\link{Kinhom}}. If 
#' \code{renormalise=T}, the estimated intensity \eqn{lambda} is multiplied by
#' \eqn{c^(normpower/2)}, where \eqn{c = area(W)/sum[i] (1/lambda(x[i]))}. This way,
#' renormalisation has about the same effect as in \code{\link{Kinhom}}.
#' 
#' All edge corrections are implemented as approximations. The translational edge 
#' correction suffers from a small intrinsic bias in some cases of the locally 
#' scaled point processes, depending on intensity and window shape. It is 
#' recommended to use Ripley's isotropic edge correction, if possible.
#' 
#' For more details, see \code{\link{Kscaled}}
#' @export
#' @seealso \code{\link{Kest}}, \code{\link{Kinhom}}, \code{\link{Kscaled}}


Klocscaled <- 
  function (X, lambda=NULL, ..., r = NULL, breaks = NULL, 
            correction=c("border", "isotropic", "translate"),
            renormalise=FALSE,
            normpower=1,
            sigma=NULL, varcov=NULL,
            r.max = 10)
{
    verifyclass(X, "ppp")
    rfixed <- !missing(r) || !missing(breaks)

    # determine basic parameters
    W <- X$window
    npts <- npoints(X)
    area <- area.owin(W)

    # match corrections
    correction.given <- !missing(correction) && !is.null(correction)
    correction <- pickoption("correction", correction,
                             c(none="none",
                               border="border",
                              # "bord.modif"="bord.modif",
                               isotropic="isotropic",
                               Ripley="isotropic",
                               trans="translate",
                               translate="translate",
                               translation="translate",
                               best="best"),
                             multi=TRUE)

    best.wanted <- ("best" %in% correction)
    correction <- implemented.for.K(correction, W$type, correction.given)

   ###########################################################
    # DETERMINE WEIGHTS AND VALIDATE
    #

    if(missing(lambda)) {
      # No intensity data provided
      # Estimate density by leave-one-out kernel smoothing
      lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                        at="points", leaveoneout=TRUE)
      lambda <- as.numeric(lambda)
    } else {
        # lambda values provided
      if(is.im(lambda)) 
        lambda <- safelookup(lambda, X)
      else if(is.function(lambda)) 
        lambda <- lambda(X$x, X$y)
      else if(is.numeric(lambda) && is.vector(as.numeric(lambda)))
        check.nvector(lambda, npts)
      else stop(paste(sQuote("lambda"),
                      "should be a vector, a pixel image, or a function"))
    }
  
    ################ CHANGE (U, 01.08.2013)   ####################
    # renormalise. Here we only need half the power ;-)
      if(renormalise) {
        check.1.real(normpower)
        stopifnot ((1 <= normpower) & (normpower <= 2)) #((0.5 <= normpower) & (normpower <= 1))
        renorm.factor <- (area/sum(1 / lambda))^(normpower / 2)
        lambda <- lambda / renorm.factor
      } 
    
 ################ CHANGE (U, 01.08.2013)   ####################
  # recommended range of r values
  # use max lambda instead of npts/area 
  # calculate recommended arguments r differently
  # 
    area <- area.owin(W)
     
    maxrescale <- sqrt(max(lambda))
    minrescale <- sqrt(min(lambda))
    
    rmaxdefault <- rmax.rule("K", W) 
    if(!is.null(r)) suggested.absr <- r/maxrescale 
      else suggested.absr <-NULL
    breaks <- handle.r.b.args(suggested.absr, breaks, W, 
                    rmaxdefault=rmaxdefault)
    absrmax <- breaks$max
    breaks$val <- breaks$val * maxrescale
    breaks$r <- breaks$r * maxrescale
    r <- breaks$r
    if (max(r)> r.max) 
    {r <- seq(0, r.max, length.out=513); breaks$val <-c(-r[2],r)}
    
    
    alim <- c(0, min(absrmax*maxrescale, rmaxdefault* maxrescale, r.max))
    rlim <- diameter(W)/2 * maxrescale # minrescale # this is picky, really!
  
     
# no efficient border correction engine here...     
        
    # this will be the output data frame
    K <- data.frame(r=r, theo= pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson %s")
  
   # Khatlab <- substitute(hat(K)[0]^(s)*(r), NULL)
    Ktheolab <- substitute(K[0]^(s)*(r), NULL)
    Kname <- "K[0]^{(s)}"
    Krname <- "K[0]^{(s)}(r)"
    Khatname <- "widehat(K)[0]^(s)"
    K <- fv(K, "r", Ktheolab,
            "theo", , alim, c("r", Krname), desc, fname=Kname)
 #   K <- tweak.fv.entry(K, "theo", new.labl="paste(K[0]^{(s)}*(r),', ',scriptstyle(CSR))") 
    K <- tweak.fv.entry(K, "theo", new.labl="paste(K[0]*(r),', ',scriptstyle(CSR))") 
    ############### END CHANGE ##################################   
    
    # identify all close pairs
    rmax <- max(r)
 ################ CHANGE (U, 16.08.2013)   ####################
    close <- lsclosepairs(X, lambda, rmax)
    dIJ <- close$d
    eudIJ <- close$eud
 ################ END CHANGE ##################################       
    # compute weights for these pairs
    I <- close$i
    J <- close$j
    XI <- X[I]

   if(any(correction == "none")) {
    # uncorrected! For demonstration purposes only!
    wh <- whist(dIJ, breaks$val)  # no weights
    Kun <- cumsum(wh)/npts
    K <- bind.fv(K, data.frame(un=Kun), "%s[un](r)",
                 "uncorrected estimate of %s",
                 "un")
  }
  
 if(any(correction == "border")) {
  # border method
  # Compute SCALED distances to boundary
    b <- bdist.points(X) * sqrt(lambda)
    I <- close$i
    bI <- b[I]
  # apply reduced sample algorithm
    RS <- Kount(dIJ, bI, b, breaks)
  # ################ CHANGE (U, 01.08.2013)   ####################
    Kb <- RS$numerator/RS$denom.count
    Kb[r > rlim] <- NA 
    K <- bind.fv(K, data.frame(border=Kb),"%s*(r)", # "%s[bord](r)",
                   "border-corrected estimate of %s",
                   "border")
    K <- rebadge.fv(K, Ktheolab, Khatname)
    K <- tweak.fv.entry(K, "border", new.labl="paste(%s*(r),', ',scriptstyle(bord))") 
   
   ################ END CHANGE ##################################       
  }

  if(any(correction == "translate")) {
    # translation correction
    XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    wh <- whist(dIJ, breaks$val, edgewt)
    Ktrans <- cumsum(wh)/npts 
  # ################ CHANGE (U, 01.08.2013)   ####################
  #  h <- diameter(W)/2
  #  Ktrans[r >= h] <- NA
    Ktrans[r > rlim] <- NA
   ################ END CHANGE ##################################       
    K <- bind.fv(K, data.frame(trans=Ktrans), "%s[trans](r)",
                 "translation-corrected estimate of %s",
                 "trans")
    K <- rebadge.fv(K, Ktheolab, Khatname)
    K <- tweak.fv.entry(K, "trans", new.labl="paste(%s*(r),', ',scriptstyle(trans))") 
   
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction (using UN-SCALED distances)
    edgewt <- edge.Ripley(XI, matrix(eudIJ, ncol=1))
    wh <- whist(dIJ, breaks$val, edgewt)
    Kiso <- cumsum(wh)/npts 
  ################ CHANGE (U, 01.08.2013)   ####################
  #  h <- diameter(W)/2
  #  Kiso[r >= h] <- NA
    Kiso[r > rlim] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "%s[iso](r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso")
     K <- rebadge.fv(K, Ktheolab, Khatname)
    K <- tweak.fv.entry(K, "iso", new.labl="paste(%s*(r),', ',scriptstyle(iso))") 
  }
  # ################ CHANGE (U, 01.08.2013)   ####################
  K$theo[r > rlim] <- NA    
    # compute edge corrected estimates
 
  # default plot will display all edge corrections
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  fvnames(K, ".") <- nama[!(nama %in% c("r", "rip", "ls"))]
  #
  unitname(K) <- c("normalised unit", "normalised units")
  return(K)
}

