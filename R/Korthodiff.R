# scaled version of Ohser-Stoyan directional K-function, 31.08.12 / 25.08.13 uh

# now with lambda / scalefaktor marked point processes



#' Contrast of orthogonal directional K-functions
#'
#' Estimates the difference of vertical and horizontal directional \eqn{K}-function of a point process.
#' 
#'
#' @param X a point pattern, object of class \code{"ppp"}.
#' @param type character, the type of second-stationarity assumed: 
#'   \itemize{
#'        \item \code{"s"} locally rescaled
#'        \item \code{"hs"} homogeneous, but evaluated as scaled \eqn{K}-function
#'   }
#'   Only the first match is used, and only local scaling is supported so far.
#' @param dphi half the angle of the \eqn{K}-functions directions
#' @param r optional: vector of argument values \eqn{r} at which \eqn{K(r)} should
#'   be evaluated.
#' @param correction a character vector giving the edge correction type, so far only
#'    \code{"isotropic"} is supported.
#' @param normpower an integer between 0 and 2. If \code{normpower} > 0, the 
#'  intensity is normalized, see the Details.
#' @param ... optional arguments passed to \code{\link{as.sostpp}}   
#' @param rescaleangle logical, indicating whether to rescale the directional 
#'   \eqn{K}-functions as to represent the full circle
#' @param max.ls.r upper limit for argument \eqn{r}. The default 3 is quite large - note
#'  that the argument of the locally rescaled \eqn{K}-function corresponds to
#'  \deqn{r/\sqrt{\lambda}}{r / sqrt(lambda)} in the non scaled case.
#'
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @seealso \code{\link{Kest}}, \code{\link{Kinhom}}, \code{\link{Kscaled}}
#' @examples
#' # The bronzefilter pattern is quite isotropic. Compare two halves:
#' data(bronzefilter)
#' X <- rescaled(bronzefilter)                 
#' W1 <- owin(c(0,9), c(0,7))
#' W2 <- owin(c(9,18), c(0,7))
#' rr <- seq(0, 1.3, .01)
#' plot(DeltaKdir (X[W1], r=rr), ylim=c(-2,1), col= c("red", "black"), main="")
#' lines(rr, DeltaKdir (X[W2],r = rr)$iso, col="blue") 
#'
#' # the backtransformed pattern shows anisotropy
#' Y <- backtransformed(retransformed(bronzefilter, trafo = "gradx"))                 
#' lines(rr, DeltaKdir (Y[W1], r = rr)$iso,  col= c("red"), lty = "dotted")
#' lines(rr, DeltaKdir (Y[W2], r = rr)$iso, col="blue", lty = "dotted" ) 

estDeltaKdir <- function (X, 
                       type = c("s", "t", "w", "h", "hs"),      
                       dphi = pi/4, # half angle!!! in contrast to the above
                       r = NULL,
                       correction = "isotropic",# c("border", "isotropic", "Ripley", "translate"),
                       normpower = 0,
                       ...,
                       rescaleangle = TRUE, # rescale to full angle
                       max.ls.r = 3.0) 
{  
  npts <- npoints(X)
  W <- X$window
  area <- area.owin(W)
  stopifnot(npts > 1)
  
  if (missing(type))   
  {
    sostype <- currenttype(X)
    if (length(sostype) == 0) 
      {
         X <- ashomogeneous(X)
         sostype <- "h"
      }
  }
  else {
    sostype <- type[1]
    if (!has.type (X, sostype)) X <- as.sostpp(X, type = sostype, ...)
  }
    
  marx <- X$typemarks
  if (normpower != 0) {
    stopifnot ((1 <= normpower) & (normpower <= 2)) 
    if (!is.null(marx$lambda)){  
      renorm.factor <-  (sum(1 / marx$lambda) / (area.owin(X)))^(normpower / 2) 
      marx$lambda <- marx$lambda * renorm.factor
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
  
  if (scaling) marx$lambda <- rep(1, npts) # refers to unit rate template
  
  if (sostype == "t") 
  {
    X <- backtransformed (X)
    # make uniform lambdas
    marx$lambda <- npts / area
    # approximate window, if not a rectangle
    W <- X$window
    if (!is.rectangle(W))
    {
      warning("window of backtransformed pattern approximated by its enclosing rectangle")
      W <- owin(W$xrange, W$yrange)
      X$window <- W
    }
  }
  
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
    if(sostype == "w") lamax <- max(marx$lambda) else lamax <- npts / area
    rmaxdefault <- rmax.rule("K", W, lamax)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max
    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
  }
  
  ### formalities ---------------------------------------------------------------------------------
  # name of the function and its estimates
  typename = sostype
  # if (sostype == "hs") typename ="s"
  # if (sostype == "h") typename ="w"
  
  Kname <- paste("Delta K[dir]^{(", typename, ")}", sep="")
  Krname <- paste("Delta K[dir]^{(", typename, ")}(r)", sep="")
  Khatname <- paste("widehat(Delta*K)[dir]^(", typename, ")", sep="")
  Ktheolab <- substitute(Delta*K[dir]^(name)*(r), list(name = typename))
  
  # start data frame with CSR
  K <- data.frame(r = r, theo = 0)
  desc <- c("distance argument r", "zero line %s")
  K <- fv(K, "r", Ktheolab,
          "theo", , alim, c("r",Krname), desc, fname=Kname)
  K <- tweak.fv.entry(K, "theo", new.labl="paste(Delta*K[dir]*(r),', ',scriptstyle(isotropy))") 
  
  ### actual calculation of  the scaled K-function  ----------------------------------------------- 
  # identify all close pairs
  rmax <- max(r)
  
  if(scaling) close <- lsclosepairs(X, rmax, invscale = marx$invscale)
  else close <- lsclosepairs(X, rmax) 
   
  # only those in the orientation interval around 0 and pi/2, symmetric     
  dir <-  with(close, atan2(dy,dx))%% (2*pi)
  dirv <- abs((dir%%pi) - pi/2)
  dirh <- pi/2 - dirv
  # change to: dir <- with(close, (atan2(dy,dx) - phi)%%(2*pi))
  innenh <- dirh <  dphi 
  innenv <- dirv <  dphi 
  close$in.hor <- innenh
  close$in.ver <- innenv
  
  close <- as.data.frame(close)[innenv|innenh,]
  
  dIJ <- close$d
  eudIJ <- close$eud
  # compute weights for these pairs
  I <- close$i
  J <- close$j
  XI <- X[I]
  XJ <- X[J]
  
  if (weighted) wIJ <- 1 / (marx$lambda[J] * marx$lambda[I])
  else  wIJ <- 1 / marx$lambda[J] * area / (npts - 1)

  
  edgewts <- edgecross.Ripley(XI, r = matrix(eudIJ, ncol = 1), dphi = dphi)
  if(!is.data.frame(edgewts)) stop(paste("uhadada something wrong with",X$n,"points"))
  edgewt <- close$in.hor * edgewts$w.horiz - close$in.ver * edgewts$w.vert
  totalwt <- edgewt * wIJ
  wh <- whist(dIJ, breaks$val, totalwt)
  Kiso <- cumsum(wh) / area
  if (rescaleangle) Kiso <- Kiso * pi / dphi /2
  K <- bind.fv(K, data.frame(iso=Kiso), "%s[iso](r)",
               "cross Ripley isotropic correction estimate of %s", "iso")
  K <- rebadge.fv(K, Ktheolab, Khatname)
  K <- tweak.fv.entry(K, "iso", new.labl="paste(%s*(r),', ',scriptstyle(iso))") 
      
  
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  nama <- nama[nama != "r"] #!(nama %in% c("r", "rip", "ls"))]
  fvnames(K, ".") <- nama
  attr(K, "sostype") <- sostype
  if (typename == "s") unitname (K) = NULL else unitname(K) <- unitname(X)
  
  return(K) 
}


