# second-order stationarity typed point patterns

#' Convert an R-object  to class sostyppp
#'
#' Converts an R-object into a second-order
#' stationarity typed point pattern of class \code{"sostyppp"}.
#'
#' @export
#' @param x an R-object, currently of class \code{"ppp"} of \code{"sostyppp"}
#' @param type character giving the type of second order stationarity. One of
#'    "w", "t", "s", "h", "hs" (or "none").
#' @param ... further arguments
#' @return an object of class \code{"sostyppp"}, but with no type information.
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

as.sostyppp <- function(x, type = "h", ...)
{
  UseMethod("as.sostyppp")
}


#' Convert object of class ppp to class sostyppp
#'
#' Converts a spatstat-point pattern of class \code{"ppp"} into a second-order
#' stationarity typed point pattern of class \code{"sostyppp"}, and assigns type of
#' second-order stationarity.
#'
#' @S3method as.sostyppp ppp
#' @method as.sostyppp ppp
# @export
#' @param x an object of class \code{"ppp"}
#' @param type character giving the type of second order stationarity. One of
#'    "w", "t", "s", "h", "hs" (or "none").
#' @param ... arguments passed to \code{\link{reweighted}}, \code{\link{retransformed}},
#' \code{\link{rescaled}}. Ignored if \code{type} = \code{"h"} or \code{type} = \code{"hs"}.
#' @details Type \code{"none"} is only for internal use, and may result in errors when
#' trying to estimate summary statistics or alike.
#' @return an object of class \code{"sostyppp"}, but with no type information.
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' # print an ordinary ppp point pattern
#' print(bronzefilter)
#' brofi <- as.sostyppp (bronzefilter)
#' print(brofi)

as.sostyppp.ppp <- function(x, type = "h", ...)
{
  X <- x
  firstclass(X) <-  "sostyppp"
  X$sostinfo <- list()
  if (type == "none") { X$sostinfo$tmarks <- NULL;  X$sostype <- .settype("none", NULL) }
  else { if (type == "w") X <- reweighted(X, ...)
    else { if (type == "t") X <- retransformed(X, ...)
      else { if (type == "s") X <- rescaled(X, ...)
        else { if (type == "h") X <- homogeneous(X, "h")
          else { if (type == "hs") X <- homogeneous(X, "hs")
            else { warning("unknown 2ndorder stationarity type") }}}}}}
  return(X)
}

#' Convert object of class sostyppp to class sostyppp
#'
#' Wrapper function for \code{\link{reweighted}}, \code{\link{retransformed}},
#' \code{\link{rescaled}} and \code{\link{homogeneous}}, assigns type of
#' second-order stationarity.
#'
#' @S3method as.sostyppp sostyppp
#' @method as.sostyppp sostyppp
# @export
#' @param x an object of class \code{"sostyppp"}
#' @param type character giving the type of second order stationarity. One of
#'    "w", "t", "s", "h", or "hs".
#' @param ... arguments passed to \code{\link{reweighted}}, \code{\link{retransformed}},
#' \code{\link{rescaled}}. Ignored if \code{type} = \code{"h"} or \code{type} = \code{"hs"}.
#' @return an object of class \code{"sostyppp"} with corresponding type information.
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' brofi <- as.sostyppp (bronzefilter)
#' bronze <- as.sostyppp (brofi, "t", "gradx")

as.sostyppp.sostyppp <- function(x, type = "h", ...)
{
  if (type == "none") {
    x$sostinfo$tmarks <- NULL
    x$sostype <- .settype("none", NULL)
  }
  else if (type == "w") x <- reweighted(x, ...)
    else { if (type == "t") x <- retransformed(x, ...)
      else { if (type == "s") x <- rescaled(x, ...)
        else { if (type == "h") x <- homogeneous(x, "h")
          else { if (type == "hs") x <- homogeneous(x, "hs")
            else { warning("unknown 2ndorder stationarity type") }}}}}
  return(x)
}


#' Make point pattern reweighted 2nd-order stationary
#'
#' Add type information for reweighted second-order stationary point processes.
#'
#' @param X the original point pattern, an object of class \code{"ppp"} or \code{"sostyppp"}
#' @param intensity optional. The estimated intensity function, can be either
#' \itemize{
#'   \item a single number,
#'   \item a vector of intensity values at the points of \code{X},
#'   \item a \code{function(x,y)} that returns the intensity and is valid at each point of \code{X}
#'   \item a pixel image of class \code{"\link{im}"}
#'   }
#' @param \ldots optional extra parameters. If \code{intensity} is given as a function, \ldots may contain extra parameters,
#'   if \code{intensity} is empty, these are parameters passed to \code{spatstat::\link[spatstat]{density.ppp}}.
#' @return A reweighted s.o.s. typed point pattern (object of class \code{{"sostyppp"}}),
#' having typemarks  with  an element \code{intensity} (estimated intensity)
#' @export
#' @details
#'   If  \code{intensity} is missing, the function will check if the pattern \code{X}
#'   has previously been typed as rescaled s.o.s. (by function \code{\link{rescaled}}). In that case,
#'   the intensity is calculated as the square inverse scale factor.
#'
#'   If  \code{intensity} is empty, and \code{X} has no inverse scale factor marks, the intensity
#'   is estimated using \code{spatstat::\link[spatstat]{density.ppp}}, and extra arguments \ldots
#'   are passed to \code{density.ppp}. By contrast to \bold{spatstat},
#'   the {"leaveoneout"}-method is not chosen by default.
#'
#'   If  \code{intensity} is given as a function, extra parameters may be passed as \ldots.
#'
#'   If \code{intensity} is given as a single number, the point pattern gets marked with
#'   constant intensity, which effectively means that it will be analysed as a homogeneous
#'   point process with known intensity.
#'
#' @seealso related functions: \code{\link{rescaled}}, \code{\link{retransformed}}, \code{\link{homogeneous}}
#' @seealso \code{\link{sostyppp}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

reweighted <- function (X, intensity = NULL, ...)#, normpower = 0)
{
  if(!is.sostyppp(X)) X <- as.sostyppp.ppp(X, "none")
  npts <- npoints(X)
  marx <- X$sostinfo$tmarks
  if (npts > 0){
    if (is.null(intensity) && !is.null(marx$invscale)) intensity <- marx$invscale^2
    la <- getIntensity(X, intensity, ...)#,  normpower = normpower)
    if (!is.null (marx$intens)) marx$intens <- la
    else marx <- as.data.frame(cbind(marx, intens = la))
    # now marx is really a data frame
    X$sostinfo$tmarks <-  marx
  }
 # X$sostinfo$intensity <- intensity
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
#' @param intensity optional extra parameters, see details.
#' The possible format of this corresponds to that of parameter \code{invscale}.
#' @param \ldots arguments passed to \code{invscale} or \code{intensity}, if these are functions, or to
#'   \code{spatstat::\link[spatstat]{density.ppp}} if both \code{invscale} and \code{intensity} are missing.
# @param normpower an integer between 0 and 2. If \code{normpower} > 0 and \code{invscale} is missing,
#   the intensity \code{intensity} is normalized,
#   see the Details.
#' @return A rescaled s.o.s. typed point pattern (object of class \code{"\link{sostyppp}"}), having typemarks
#'   with   an element  \code{invscale} (square root of estimated intensity)
#' @details
#'   If \code{invscale} is empty, the square root of the intensity \code{itensity} is used instead.
#'   If neither \code{invscale} nor \code{intensity} are given, the function checks if the pattern \code{X}
#'   has previously been marked with the intensity (by function \code{\link{reweighted}}).
#'   If both \code{invscale} and \code{intensity} are empty and there are no intensity marks, the intensity
#'   is estimated using \code{spatstat::\link[spatstat]{density.ppp}}, and extra arguments \ldots
#'   are passed to \code{density.ppp}.
#'
#'   If \code{invscale} (or \code{intensity}) is given as a function, extra parameters may be passed as \ldots.
#'
#'   If \code{invscale} is given as a single number, the point pattern gets marked with
#'   constant scale factor, which effectively means that it will be analysed as a homogeneous
#'   point process.
#' @seealso related functions: \code{\link{reweighted}}, \code{\link{retransformed}}, \code{\link{homogeneous}}
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#'@references
#'Hahn, U. and Jensen, E. B. V. (2013)
#'  Inhomogeneous spatial point processes with hidden second-order stationarity.
#'  \emph{CSGB preprint} 2013-7.
#'  \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#'

rescaled <- function (X,  invscale = NULL, intensity = NULL, ...)
{
  if(!is.sostyppp(X)) X <- as.sostyppp.ppp(X, "none")

  npts <- npoints(X)
  marx <- X$sostinfo$tmarks

  if (npts > 0){
  if (is.null(invscale))  {
    if (is.null(intensity)){
       if (!is.null(marx$intens)) la <- marx$intens
       else la <- getIntensity(X, NULL, ...)
    }
    else la <- getIntensity(X, intensity, ...)
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
   X$sostinfo$tmarks <-  marx
  }
#  X$sostinfo$invscale <- invscale
#  X$sostinfo$intensity <- intensity
  X$sostype <- .settype("s", X$sostype)
  return(X)
}

#' Make point pattern retransformed 2nd-order stationary
#'
#' Add mark information for retransformed second-order stationary point processes.
#'
#' @param X the original point pattern
#' @param backtrafo the coordinate transformation, can be either
#' \itemize{
#'   \item a \code{function(x,y)} that returns a data frame or list with elements \code{x} and \code{y}
#'    and is valid at each point of \code{X}
#'   \item a string with value \code{"gradx"} or \code{"grady"}.
#'   }
#' @param intensity optional, intensity. Only used if \code{backtrafo="gradx")} or \code{backtrafo="grady")}.
#' can be given as  \itemize{
#'   \item a single number,
#'   \item a vector of intensity values at the points of \code{X},
#'   \item a \code{function(x,y)} that returns the intensity and is valid at each point of \code{X}
#'   \item a pixel image of class \code{"\link{im}"}.
#'   }
#' @param trafo optional; a \code{function(x,y)} that returns a data frame or list with elements \code{x} and \code{y}
#'    and is valid at each point of \code{X}. The inverse of \code{backtrafo},
#' is used by function \code{"\link{backtransformed}"} to transform the window if it is a binary mask.
#' @param ... optional extra parameters for \code{backtrafo}, if this is given as a \code{function}.
#' @return A retransformed s.o.s. typed point pattern (object of class \code{"\link{sostyppp}"})
#' with information about the backtransform.
#' @seealso related functions: \code{\link{reweighted}}, \code{\link{rescaled}}, \code{\link{homogeneous}}
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @details
#'   If \code{backtrafo = "gradx"} or \code{backtrafo = "grady"}, the backtransformation is obtained
#'   from the  data by inverse cdf transform, affecting only the \eqn{x}- or \eqn{y}-coordinate,
#'   respectively. For these gradient transformations, it is assumed that the
#'   transformation is one-to-one, in the sense that the observation window is mapped to itself.
#'   If additionally the intensity \code{intensity} is given, the gradient back transform is
#'   calculated from \code{intensity} as described in Hahn & Jensen (2013).
#'
#'   The \code{tinfo} element in the resulting \code{sostyppp}-object contains two functions
#'   \code{transform} and \code{backtransform}. They are used when the pattern is backtransformed to
#'   homogeneity by function with \code{"\link{backtransformed}"}.
#'   In the case of gradient transformations, these functions are obtained by linear transformation
#'   with (\code{\link{approxfun}}). The reverse of \code{backtransform} is needed when the window is
#'   a binary mask, see \code{spatstat::\link[spatstat]{owin}}, and has to be given as argument {\code{trafo}}.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#'@references
#'Hahn, U. and Jensen, E. B. V. (2013)
#'  Inhomogeneous spatial point processes with hidden second-order stationarity.
#'  \emph{CSGB preprint} 2013-7.
#'  \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#
# look at akima, function interp

retransformed <- function (X, backtrafo = identxy, intensity = NULL, trafo = NULL, ...)
{
  if(!is.sostyppp(X)) X <- as.sostyppp.ppp(X, "none")
  npts <- npoints(X)
  isGradTrafo <- (backtrafo %in% c("gradx", "grady"))
  gradient <- ifelse(isGradTrafo, backtrafo, NULL)
  if (npts > 0){
    if (is.function(backtrafo)){
      if (is.mask(X$window)) if (!is.function(trafo))
         warning("retransformed point pattern with binary mask, but function trafo not given")
    }

# gradient transform

    else {
      if (isGradTrafo) {
        if (is.null(intensity)){  # just an ordinary inverse cdf transformation
          if (backtrafo == "gradx") {
            xlims <- X$window$xrange
            oo <- rank(X$x, ties.method = "first")
            x0 <- (xlims[1] + diff(xlims) * (0.5 + (1:npts)) / (npts+1))[oo]
            y0 <- X$y
          }
          else {
            ylims <- X$window$yrange
            oo <- rank(X$y, ties.method = "first")
            y0 <- (ylims[1] + diff(ylims) * (0.5 + (1:npts)) / (npts+1))[oo]
            x0 <- X$x
          }
        }
        else {
          la <- getIntensity(X, intensity, ...)
          if(backtrafo == "gradx") {
            z <- X$x; zrange <- X$window$xrange
          }
          else {
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
          if (backtrafo == "gradx") {
            x0 <- z
            y0 <- X$y
          } else {
            y0 <- z
            x0 <- X$x
          }
        }

        if (backtrafo == "gradx")
        {
          xw <- X$window$xrange
          xx0 <- c(xw[1], x0, xw[2])
          xx <- c(xw[1], X$x, xw[2])
          trafo <- function(x, y = NULL) {
            if (is.null(y)) return(data.frame(x = approxfun(xx0, xx, ties = mean)(x$x), y = x$y))
            else return (data.frame(x = approxfun(xx0, xx, ties = mean)(x), y = y))
          }
          backtrafo  <- function(x, y = NULL) {
            if (is.null(y)) return(data.frame(x = approxfun(xx, xx0, ties = mean)(x$x), y = x$y))
            else return (data.frame(x = approxfun(xx, xx0, ties = mean)(x), y = y))
          }
        }
        else
        {
          yw <- X$window$yrange
          yy0 <- c(yw[1], y0, yw[2])
          yy <- c(yw[1], X$y, yw[2])
          trafo  <- function(x, y = NULL) {
            if (is.null(y)) return(data.frame(x = x$x, y = approxfun(yy0, yy)(x$y)))
            else return (data.frame(x = x, y = approxfun(yy0, yy)(y)))
          }
          backtrafo  <- function(x, y = NULL) {
            if (is.null(y)) return(data.frame(x = x$x, y = approxfun(yy, yy0)(x$y)))
            else return(data.frame(x = x, y = approxfun(yy, yy0)(y)))
          }
        }
      }
      else stop(paste(sQuote("trafo"),
                      "should be a function, or one of the character strings",
                      dQuote("gradx"),"or", dQuote("grady")))
  }
  }
  X$sostinfo$backtransform <- backtrafo
  X$sostinfo$transform <- trafo
  # X$sostinfo$isGradTrafo <- isGradTrafo
  X$sostinfo$gradient <- gradient
  X$sostype <- .settype("t", X$sostype)
#  X$sostinfo$intensity <- intensity
  return(X)
}

#' Mark point pattern as homogeneous
#'
#' Add type information for second-order stationary point processes.
#'
#' @param X the original point pattern
#' @param type optional character. If set to "hs", the patterns is preferentially
#'   evaluated as "rescaled", otherwise as "reweighted" (default).
#' @param intensity optional, constant intensity.
#' @return A homogeneous s.o.s. typed point pattern (object of class \code{"\link{sostyppp}"}), having typemarks
#'   with   elements \code{intensity} (estimated intensity) and
#'   \code{invscale} (square root of intensity).
#' @seealso related functions: \code{\link{rescaled}}, \code{\link{retransformed}}, \code{\link{reweighted}}
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

homogeneous <- function (X,  type="h", intensity = NULL)
{
  if(!is.sostyppp(X)) X <- as.sostyppp.ppp(X, "none")
  npts <- npoints(X)
  marx <- X$sostinfo$tmarks
  if (npts > 0){
    if (is.null(intensity)) la <- npts / area.owin(X$window) else la <- intensity[1]
    if (!is.null(marx$intens)) marx$intens <- la
    else marx <-  as.data.frame(cbind(marx, intens = rep(la, npts)))
    if (!is.null(marx$invscale)) marx$invscale <- sqrt(la)
    else marx <-  cbind(marx, invscale = sqrt(la))
    # make sure last type is set to given argument:
    X$sostinfo$tmarks <-  marx
  }
#  X$sostinfo$intensity <- intensity
  X$sostype  <- .settype(type, .settype("h", .settype("hs", NULL)))
  return(X)
}


#' Obtain backtransformed point pattern
#'
#' Backtransform a pattern that is tagged as  retransformed second-order stationary.
#'
#' @param X a retransformed s.o.s. type tagged point pattern, object of class \code{{sostyppp}}.
#' @return a homogeneous s.o.s. type tagged point pattern of class \code{{sostyppp}}
#'  with points in \code{x0} and \code{y0}.
#'  The type of the result is set to both  homogeneous and scaled-homogeneous,
#'  and is typemarked with constant intensity and constant inverse scale factor.
#'
#'  The window of the resulting pattern is polygonal if the original was a rectangle or polygonal;
#'  it is a binary mask if the original window was a binary mask.
#' @details As retransformed s.o.s. typed point pattern, \code{X} has an element \code{tinfo},
#' a list containing information about the transformation. If this list contains
#' an element \code{backtransform}, this element is taken to be a function and used
#' to backtransform the window. If the original transformation was a gradient transformation,
#' the coordinates of the window are interpolated. Currently, no saftey precautions if
#' \code{backtransform} returns rubbish.
#'
# If no element \code{invtrafo} is present in \code{X$extra}, the original window is retained, which also might result in rubbish.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' bronzetra <- retransformed(bronzefilter, "gradx")
#' bronzetemplate <- backtransformed(bronzetra)
#' plot(bronzetemplate, use.marks = FALSE)
backtransformed <- function(X)
{
  stopifnot(is.sostyppp(X))
  stopifnot(hasType(X, type= "t"))
  if (is.null(X$sostinfo$backtransform)) stop ("no backtransform given")
  preserveRectangle <- is.rectangle(X$window) & !is.null(X$sostinfo$gradient)
  if (is.function(X$sostinfo$transform))
    Y <- coordTransform(as.ppp(X),
                        trafoxy = X$sostinfo$backtransform,
                        invtrafoxy = X$sostinfo$transform,
                        subdivideBorder = !preserveRectangle)
  else Y <- coordTransform(as.ppp(X),
                           trafoxy = X$sostinfo$backtransform,
                           subdivideBorder = !preserveRectangle)
  return(homogeneous(Y))
}



# @param X point pattern, of class ppp
# @param intensity optional intensity, number, vector, function or image
# @param ... extra parameters for estimation of intensity
#' @rdname sostatpp-internal
#' @export
#' @keywords internal

getIntensity <- function(X, intensity = NULL, ...)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if(is.null(intensity)) {
    # No intensity data provided
    # use a kernel estimate
    # in spatstat, leaveoneout = TRUE is set by default. This may lead to strange
    # result in sparse regions, in particular when the density is "renormalized"
    # afterwards. Therefore this default has been removed here.
    intensity <- density(X, ..., at="points")
    intensity <- as.numeric(intensity)
  } else {
    # intensity values provided
    if(is.im(intensity))
      intensity <- safelookup(intensity, X)
    else if(is.function(intensity))
      intensity <- intensity(X$x, X$y, ...)
    else if(is.numeric(intensity) && is.vector(as.numeric(intensity))){
      if(length(intensity)==1) intensity <- rep(intensity, npts)
      else check.nvector(intensity, npts)
    }
    else stop(paste(sQuote("intensity"),
                    "should be a vector, a pixel image, or a function"))
  }
  return(intensity)
}

# renormalize intensity as in spatstat
# @param X point pattern, of class ppp or sostyppp
# @param intensity optional intensity vector
# @param normpower renormalizing power.
# @details if intensity is not given, and X is sos retransformed or reweighted,
# the typemarks are adjusted accordingly
#' @export
#' @rdname sostatpp-internal
#' @keywords internal

normalizedIntensity <- function(X, intensity = NULL, normpower = 2)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if(is.null(intensity)) {
    if(!is.sostyppp(X)) stop("no intensity given")
    # first priority: current type
    if (currentType(X) %in% c("w","s")) {
        if (currentType(X) == "w") intensity <- X$sostinfo$tmarks$intens
        else intensity <- X$sostinfo$tmarks$invscale^2
      }
    else if (hasType(X, "w")) intensity <- X$sostinfo$tmarks$intens
    else if (hasType(X, "s")) intensity <- X$sostinfo$tmarks$invscale^2
    else stop("error: no intensity information in point pattern")
  }
  else {
    if(length(intensity) == 1) intensity <- rep(intensity, npts)
    else if(length(intensity) != npts)
      stop(paste(sQuote("intensity"),
        "should be a single number or a vector of values for each point in",
        sQuote("X")))
  }
  if (normpower != 0) {
    stopifnot ((1 <= normpower) & (normpower <= 2))
    renorm.factor <-  (sum(1 / intensity) / (area.owin(X)))^(normpower / 2)
    intensity <- intensity * renorm.factor
  }
  return(intensity)
}