# second-order stationarity typed point patterns

#' Convert an R-object  to class sostpp
#' 
#' Converts an R-object into a second-order 
#' stationarity typed point pattern of class \code{"sostpp"}. 
#'
#' @export
#' @param x an R-object, currently of class \code{"ppp"} of \code{"sostpp"}
#' @param type character giving the type of second order stationarity. One of
#'    "w", "t", "s", "h", "hs" (or "none").
#' @param further arguments    
#' @return an object of class \code{"sostpp"}, but with no type information.
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

as.sostpp <- function(x, type = "h", ...) 
{
  UseMethod("as.sostpp")
}


#' Convert object of class ppp to class sostpp
#' 
#' Converts a spatstat-point pattern of class \code{"ppp"} into a second-order 
#' stationarity typed point pattern of class \code{"sostpp"}, and assigns type of
#' second-order stationarity.
#'
#' @S3method as.sostpp ppp
#' @export
#' @param x an object of class \code{"ppp"}
#' @param type character giving the type of second order stationarity. One of
#'    "w", "t", "s", "h", "hs" (or "none").
#' @param ... arguments passed to \code{\link{reweighted}}, \code{\link{retransformed}},
#' \code{\link{rescaled}}. Ignored if \code{type} = \code{"h"} or \code{type} = \code{"hs"}.
#' @details Type \code{"none"} is only for internal use, and may result in errors when
#' trying to estimate summary statistics or alike.
#' @return an object of class \code{"sostpp"}, but with no type information.
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' # print an ordinary ppp point pattern
#' print(bronzefilter)
#' brofi <- as.sostpp (bronzefilter)
#' print(brofi)

as.sostpp.ppp <- function(x, type = "h", ...)
{
  X <- x
  class(X) <- c("sostpp", class(x))
  
  if (type == "none") { X$typemarks <- NULL;  X$sostype <- NULL }
  else { if (type == "w") X <- reweighted(X, ...)
    else { if (type == "t") X <- retransformed(X, ...)
      else { if (type == "s") X <- rescaled(X, ...)
        else { if (type == "h") X <- ashomogeneous(X, "h")
          else { if (type == "hs") X <- ashomogeneous(X, "hs")
            else { warning("unknown 2ndorder stationarity type") }}}}}}
  return(X)
}

#' Convert object of class sostpp to class sostpp
#' 
#' Wrapper function for \code{\link{reweighted}}, \code{\link{retransformed}},
#' \code{\link{rescaled}} and \code{\link{ashomogeneous}}, assigns type of
#' second-order stationarity.
#'
#' @S3method as.sostpp sostpp
#' @export
#' @param x an object of class \code{"sostpp"}
#' @param type character giving the type of second order stationarity. One of
#'    "w", "t", "s", "h", or "hs".
#' @param ... arguments passed to \code{\link{reweighted}}, \code{\link{retransformed}},
#' \code{\link{rescaled}}. Ignored if \code{type} = \code{"h"} or \code{type} = \code{"hs"}.
#' @return an object of class \code{"sostpp"} with corresponding type information.
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' brofi <- as.sostpp (bronzefilter) 
#' bronze <- as.sostpp (brofi, "t", "gradx")

as.sostpp.sostpp <- function(x, type = "h", ...)
{
  if (type == "w") X <- reweighted(x, ...)
    else { if (type == "t") X <- retransformed(x, ...)
      else { if (type == "s") X <- rescaled(x, ...)
        else { if (type == "h") X <- ashomogeneous(x, "h")
          else { if (type == "hs") X <- ashomogeneous(x, "hs")
            else { warning("unknown 2ndorder stationarity type") }}}}}
  return(X)
}




#' Make point pattern reweighted 2nd-order stationary
#' 
#' Add type information for reweighted second-order stationary point processes.
#' 
#' @param x the original point pattern, an object of class \code{"ppp"} or \code{"sostpp"}
#' @param lambda optional. The estimated intensity function, can be either 
#' \itemize{
#'   \item a single number,
#'   \item a vector of intensity values at the points of \code{X},  
#'   \item a \code{function(x,y)} that returns the intensity and is valid at each point of \code{X} 
#'   \item a pixel image of class \code{"\link{im}"}
#'   }
#' @param \ldots optional extra parameters. If \code{lambda} is given as a function, \ldots may contain extra parameters,
#'   if \code{lambda} is empty, these are parameters passed to the \bold{spatstat}-function \code{\link{density.ppp}}.
#' @return A reweighted s.o.s. typed point pattern (object of class \code{{"sostpp"}}), having typemarks 
#'   with  an element \code{lambda} (estimated intensity) 
#' @export
#' @details
#'   If  \code{lambda} is missing, the function will check if the pattern \code{X} 
#'   has previously been typed as rescaled s.o.s. (by function \code{\link{rescaled}}). In that case,
#'   the intensity is calculated as the square inverse scale factor.
#'   
#'   If  \code{lambda} is empty, and \code{X} has no inverse scale factor marks, the intensity
#'   is estimated using the \bold{spatstat}-function \code{\link{density.ppp}}, and extra arguments \ldots
#'   are passed to \code{density.ppp}. By contrast to \bold{spatstat},
#'   the {"leaveoneout"}-method is not chosen by default.
#'    
#'   If  \code{lambda} is given as a function, extra parameters may be passed as \ldots.
#'   
#'   If \code{lambda} is given as a single number, the point pattern gets marked with 
#'   constant intensity, which effectively means that it will be analysed as a homogeneous 
#'   point process with known intensity.
#'   
#' @seealso related functions: \code{\link{rescaled}}, \code{\link{retransformed}}, \code{\link{ashomogeneous}}   
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

reweighted <- function (X, lambda = NULL, ..., normpower = 0)
{
  if(!is.sostpp(X)) X <- as.sostpp.ppp(X, "none") 
  npts <- npoints(X)
  marx <- X$typemarks
  if (npts > 0){
    if (is.null(lambda) && !is.null(marx$invscale)) lambda <- marx$invscale^2
    la <- getIntensity(X, lambda, ...,  normpower = normpower)
    if (!is.null (marx$lambda)) marx$lambda <- la 
    else marx <- as.data.frame(cbind(marx, lambda = la))
    # now marx is really a data frame
    X$typemarks <-  marx
  }
  X$sostype <- .settype("w", X$sostype)
  return(X)
}

#' Make point pattern locally rescaled 2nd-order stationary
#' 
#' Add type information for locally rescaled second-order stationary point processes.
#' 
#' @param X the original point pattern
#' @param invscale optional. The \emph{inverse} of estimated scale function, can be either 
#' \itemize{
#'   \item a single number,
#'   \item a numeric vector of (inverse!) scale factor values at the points of \code{X},  
#'   \item a \code{function(x,y)} that returns the inverse scale factor and is valid at each point of \code{X} 
#'   \item a pixel image of class \code{"\link{im}"}
#'   }
#' @param lambda optional extra parameters, see details.
#' The possible format of this corresponds to that of parameter \code{iscale}.
#' @param \ldots arguments passed to \code{invscale} or \code{lambda}, if these are functions, or to 
#'   the \bold{spatstat}-function \code{\link{density.ppp}} if both \code{invscale} and \code{lambda} are missing.
# @param normpower an integer between 0 and 2. If \code{normpower} > 0 and \code{invscale} is missing, 
#   the intensity \code{lambda} is normalized,
#   see the Details.
#' @return A rescaled s.o.s. typed point pattern (object of class \code{\link{"sostpp"}}), having typemarks 
#'   with   an element  \code{invscale} (square root of estimated intensity) 
#' @details
#'   If \code{invscale} is empty, the square root of the intensity \code{lambda} is used instead.
#'   If neither \code{invscale} nor \code{lambda} are given, the function checks if the pattern \code{X} 
#'   has previously been marked with the intensity (by function \code{\link{reweighted}}).
#'   If both \code{invscale} and \code{lambda} are empty and there are no intensity marks, the intensity
#'   is estimated using the \bold{spatstat}-function \code{\link{density.ppp}}, and extra arguments \ldots
#'   are passed to \code{density.ppp}.
#'    
#'   If \code{invscale} (or \code{lambda}) is given as a function, extra parameters may be passed as \ldots.
#' 
#'   If \code{invscale} is given as a single number, the point pattern gets marked with 
#'   constant scale factor, which effectively means that it will be analysed as a homogeneous 
#'   point process.
#' @seealso related functions: \code{\link{reweighted}}, \code{\link{retransformed}}, \code{\link{ashomogeneous}}     
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}


rescaled <- function (X,  invscale = NULL, lambda = NULL, ...)
{
  if(!is.sostpp(X)) X <- as.sostpp.ppp(X) 
  
  npts <- npoints(X)
  marx <- X$typemarks
  
  if (npts > 0){
  if (is.null(invscale))  {
    if (is.null(lambda)){
       if (!is.null(marx$lambda)) la <- marx$lambda
       else la <- getIntensity(X, NULL, ...)
    }  
    else la <- getIntensity(X, lambda, ...)   
    iscale <- sqrt(la)
  } else {
      if(is.im(invscale)) 
        iscale <- safelookup(invscale, X)
      else if(is.function(invscale)) 
        iscale <- invscale(X$x, X$y, ...)
      else if(is.numeric(invscale) && is.vector(as.numeric(invscale))) {
        if (length(invscale)==1) iscale <- rep(invscale, npts) else iscale <- invscale
        check.nvector(iscale, npts)
        } 
      else stop(paste(sQuote("invscale"),
                    "should be a vector, a pixel image, or a function"))
      }         
  if (!is.null (marx$invscale)) marx$invscale <- iscale
  else marx <- as.data.frame(cbind(marx, invscale = iscale))
   X$typemarks <-  marx
  }
  X$sostype <- .settype("s", X$sostype)
  return(X)
}

#' Make point pattern retransformed 2nd-order stationary
#' 
#' Add mark information for retransformed second-order stationary point processes.
#' 
#' @param X the original point pattern
#' @param trafo the coordinate transformation, can be either 
#' \itemize{
#'   \item a \code{function(x,y)} that returns a data frame or list with elements \code{x} and \code{y}
#'    and is valid at each point of \code{X} 
#'   \item a string with value \code{"gradx"} or \code{"grady"}.
#'   }
#' @param lambda optional, intensity. Only used if \code{trafo="gradx")} or \code{trafo="grady")}. 
#' can be given as  \itemize{
#'   \item a single number,
#'   \item a vector of intensity values at the points of \code{X},  
#'   \item a \code{function(x,y)} that returns the intensity and is valid at each point of \code{X} 
#'   \item a pixel image of class \code{"\link{im}"}.
#'   }  
#' @param invtrafo optional; a \code{function(x,y)} that returns a data frame or list with elements \code{x} and \code{y}
#'    and is valid at each point of \code{X}. The inverse of \code{trafo}, 
#' is used by function \code{"\link{backtransformed}"} to transform the window. 
#' @param ... optional extra parameters for \code{trafo}, if this is given as a \code{function}.
#' @return A retransformed s.o.s. typed point pattern (object of class \code{\link{"sostpp"}}), having typemarks 
#'   with  elements \code{x0} and  \code{y0} (backtransformed points), and extra information
#'   about the backtransform. 
#' @seealso related functions: \code{\link{reweighted}}, \code{\link{rescaled}}, \code{\link{ashomogeneous}}     
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @details
#'   It is assumed that the transformation is one-to-one, in the sense that the observation window is mapped to itself.
#'   Whether this is true is \emph{not} checked. If you need to analyse transformations that are not one-to-one,
#'   this can be accomplished manually by applying the function \code{\link{coordTransform}} to the pattern \code{X},
#'   and analysing the resulting pattern as a homogeneous point pattern.
#'   
#'   If \code{trafo = "gradx"} or \code{trafo = "grady"}, the backtransformation is obtained from the 
#'   data by inverse cdf transform, under
#'   the assumption that the transformation function affects only the \eqn{x}- or \eqn{y}-coordinate, respectively.
#'   If additionally the intensity \code{lambda} is given, the gradient back transform is calculated from \code{lambda}.
#'   
#'   The function \code{invtrafo} is applied to the window when the pattern is backtransformed
#'   using \code{\link{"backtransformed"}}. If not given, linear interpolation is used, however
#'   currently only for gradient transformations (\code{trafo = "gradx"} or \code{trafo = "grady"}).
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

# look at akima, function interp

retransformed <- function (X, trafo = identxy, lambda = NULL, invtrafo = NULL, ...)
{
  if(!is.sostpp(X)) X <- as.sostpp.ppp(X, "none") 
  npts <- npoints(X)
  if (npts > 0){
    if (is.function(trafo))
    {
      xy <- trafo(X$x, X$y, ...)
      x0 <- xy$x
      y0 <- xy$y
    } 
    else if (trafo %in% c("gradx", "grady"))
    { # just an ordinary inverse cdf transformation
      if (is.null(lambda))  # just an ordinary inverse cdf transformation
      {
        if (trafo == "gradx") {
          xlims <- X$window$xrange
          oo <- rank(X$x, ties="first")
          x0 <- (xlims[1] + diff(xlims) * (0.5 + (1:npts)) / (npts+1))[oo]
          y0 <- X$y
        }
        else {
          ylims <- X$window$yrange
          oo <- rank(X$y, ties="first")
          y0 <- (ylims[1] + diff(ylims) * (0.5 + (1:npts)) / (npts+1))[oo]
          x0 <- X$x
        }
      } else {
        la <- getIntensity(X, lambda, ...)
        if(trafo == "gradx") {
          z <- X$x; zrange <- X$window$xrange
        } else {
          z <- X$y; zrange <- X$window$yrange
        }
        oo <- order(z)
        rk <- rank(z)  
        n <- length(z)  
        zz <- c(zrange[1], z[oo], zrange[2] )
        lam <- (c(la[oo], la[oo[n]]) + c(la[oo[1]], la[oo])) / 2
        intlam <- cumsum(c(0, lam * diff(zz)))
        znew <- zrange[1] + intlam / max(intlam) * diff(zrange)
        z <- znew[-c(1, n+2)][rk] 
        if (trafo == "gradx") {
          x0 <- z
          y0 <- X$y
        } else {
          y0 <- z
          x0 <- X$x
        }
      }  
    }
    else stop(paste(sQuote("trafo"),
      "should be a function, or one of the character strings",
      dQuote("gradx"),"or", dQuote("grady")))
    
    marx <- X$typemarks
    # gave a funny name for old marks in the old version using marks(X)  
    if (!is.null (marx$x0)) marx$x0 <- x0
    else marx <-  as.data.frame(cbind(marx, x0 = x0))
    if (!is.null (marx$y0)) marx$y0 <- y0
    else marx <-  cbind(marx, y0 = y0)
    X$typemarks <-  marx
   }
  # now get backtransform
  if (is.function(trafo)) 
      { itra <- trafo 
        otra <- invtrafo  }
  else{
      if (trafo == "gradx") 
      {
        xw <- X$window$xrange
        otra <- function(x, y = NULL) {
          if (is.null(y)) return(data.frame(x = approxfun(c(x0, xw), c(X$x, xw))(x$x), y = x$y))
          else return (data.frame(x = approxfun(c(x0, xw), c(X$x, xw))(x), y = y))
        }
        itra <- function(x, y = NULL) {
          if (is.null(y)) return(data.frame(x = approxfun(c(X$x, xw), c(x0, xw))(x$x), y = x$y))
          else return (data.frame(x = approxfun(c(X$x, xw), c(x0, xw))(x), y = y))
        }
      }
      else if (trafo == "grady") 
      {
        yw <- X$window$yrange
        otra <- function(x, y = NULL) {
          if (is.null(y)) return(data.frame(x = x$x, y = approxfun(c(y0, yw), c(X$y, yw))(x$y)))
          else return (data.frame(x = x, y = approxfun(c(y0, yw), c(X$y, yw))(y)))
        }
        itra <- function(x, y = NULL) {
          if (is.null(y)) return(data.frame(x = x$x, y = approxfun(c(y0, yw), c(X$y, yw))(x$y)))
          else return(data.frame(x = x, y = approxfun(c(y0, yw), c(X$y, yw))(y)))
        }
      } 
    }
  X$sostype <- .settype("t", X$sostype)
  X$extra$trafo <- trafo
  X$extra$lambda <- lambda
  X$extra$trafoo <- otra
  X$extra$invtrafo <- itra
  return(X)
}  


#' Mark point pattern as homogeneous
#' 
#' Add type information for second-order stationary point processes.
#' 
#' @param X the original point pattern
#' @param type optional character. If set to "hs", the patterns is preferentially
#'   evaluated as "rescaled", otherwise as "reweighted" (default).
#' @param lambda optional, constant lambda.   
#' @return A homogeneous s.o.s. typed point pattern (object of class \code{\link{"sostpp"}}), having typemarks 
#'   with   elements \code{lambda} (estimated intensity) and
#'   \code{invscale} (square root of lambda).
#' @seealso related functions: \code{\link{rescaled}}, \code{\link{retransformed}}, \code{\link{reweighted}}
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

ashomogeneous <- function (X,  type="h", lambda = NULL)
{
  if(!is.sostpp(X)) X <- as.sostpp.ppp(X, "none") 
  npts <- npoints(X)
  marx <- X$typemarks
  if (npts > 0){
    if (is.null(lambda)) la <- npts / area.owin(X$window) else la <- lambda[1]
    if (!is.null(marx$lambda)) marx$lambda <- la
    else marx <-  as.data.frame(cbind(marx, lambda = rep(la, npts)))
    if (!is.null(marx$invscale)) marx$invscale <- sqrt(la)
    else marx <-  cbind(marx, invscale = sqrt(la))
    # make sure last type is set to given argument:
    X$typemarks <-  marx
  }
  X$sostype  <- .settype(type, .settype("h", .settype("hs", NULL)))
  return(X)
}


#' Obtain backtransformed pattern
#' 
#' Extract the backtransformed template from a retransformed point process.
#' 
#' @param X a retransformed s.o.s. typed point pattern, object of class \code{{sostpp}}. 
#' @return a homogeneous s.o.s. typed point pattern of class \code{{sostpp}}
#'  with points in \code{x0} and \code{y0}.
#'  The type of the result is set to both  homogeneous and scaled-homogeneous,
#'  and is typemarked with constant intensity and constant inverse scale factor.
#' @details As retransformed s.o.s. typed point pattern, \code{X} has an element \code{extra}, 
#' a list containing information about the transformation. If this list contains
#' an element \code{invtrafo}, this element is taken to be a function and used 
#' to backtransform the window. If the original transformation was a gradient transformation,
#' the coordinates of the window are interpolated. Currently, no saftey precautions if 
#' \code{invtrafo} returns rubbish.
#' If no element \code{invtrafo} is present in \code{X$extra}, the original window is retained, which also might result in rubbish.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' bronzetra <- retransformed(bronzefilter, "gradx")
#' bronzetemplate <- backtransformed(bronzetra)
#' plot(bronzetemplate, use.marks = FALSE)
backtransformed <- function(X)
{
  tmarx <- X$typemarks
  stopifnot (!is.null (tmarx$x0), !is.null (tmarx$y0))
  # if inverse transform is not given, the original window is taken. This
  # can lead to problems if the transformation did not preserve the window,
  # or if the point pattern was subsampled.
  # if the trafo was gradx or grady, the result can only be approximated.
  W <- X$window
  itra <- X$extra$invtrafo
  otra <- X$extra$trafoo
  if (!is.null(itra) &!is.null(otra)) W <- coordTransform.owin(W, itra, otra)
  Y <- ashomogeneous(ppp(tmarx$x0, tmarx$y0, window = W, marks=X$marks))
  return(Y)
}



# @param X point pattern, of class ppp
# @param lambda optional intensity, number, vector, function or image
# @param ... extra parameters for estimation of intensity 
#' @rdname sostatpp-internal
#' @export
#' @keywords internal

getIntensity <- function(X, lambda = NULL, ...)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if(is.null(lambda)) {
    # No intensity data provided
    # use a kernel estimate
    # in spatstat, leaveoneout = TRUE is set by default. This may lead to strange
    # result in sparse regions, in particular when the density is "renormalized"
    # afterwards. Therefore this default has been removed here.
    lambda <- density(X, ..., at="points")
    lambda <- as.numeric(lambda)
  } else {
    # lambda values provided
    if(is.im(lambda)) 
      lambda <- safelookup(lambda, X)
    else if(is.function(lambda)) 
      lambda <- lambda(X$x, X$y, ...)
    else if(is.numeric(lambda) && is.vector(as.numeric(lambda))){
      if(length(lambda)==1) lambda <- rep(lambda, npts)
      else check.nvector(lambda, npts)
    }
    else stop(paste(sQuote("lambda"),
                    "should be a vector, a pixel image, or a function"))
  }
  return(lambda)
}

# renormalize intensity as in spatstat
# @param X point pattern, of class ppp or sostpp
# @param lambda optional intensity vector
# @param normpower renormalizing power.
# @details if lambda is not given, and X is sos retransformed or reweighted,
# the typemarks are adjusted accordingly
#' @export
#' @rdname sostatpp-internal
#' @keywords internal

normalizedIntensity <- function(X, lambda = NULL, normpower = 2)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if(is.null(lambda)) {
    if(!is.sostpp(X)) stop("no intensity given")
    # first priority: current type
    if (currenttype(X) %in% c("w","s")) {
        if (currenttype(X) == "w") lambda <- X$typemarks$lambda
        else lambda <- X$typemarks$invscale^2
      }
    else if (has.type(X, "w")) lambda <- X$typemarks$lambda
    else if (has.type(X, "s")) lambda <- X$typemarks$invscale^2
    else stop("error: no intensity information in point pattern")
  }
  else {
    if(length(lambda) == 1) lambda <- rep(lambda, npts)
    else if(length(lambda) != npts) 
      stop(paste(sQuote("lambda"),
        "should be a single number or a vector of values for each point in",
        sQuote("X")))
  }
  if (normpower != 0) {
    stopifnot ((1 <= normpower) & (normpower <= 2)) 
    renorm.factor <-  (sum(1 / lambda) / (area.owin(X)))^(normpower / 2) 
    lambda <- lambda * renorm.factor
  }
  return(lambda)
}

