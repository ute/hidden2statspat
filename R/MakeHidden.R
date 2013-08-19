# adding marks to point processes to make them hidden so stationary

#' Make point pattern reweighted 2nd-order stationary
#' 
#' Add mark information for reweighted second-order stationary point processes.
#' 
#' @param X the original point pattern
#' @param lambda optional. The estimated intensity function, can be either 
#' \itemize{
#'   \item a single number,
#'   \item a vector of intensity values at the points of \code{X},  
#'   \item a \code{function(x,y)} that returns the intensity and is valid at each point of \code{X} 
#'   \item a pixel image of class \code{"\link{im}"}
#'   }
#' @param \ldots optional extra parameters. If \code{lambda} is given as a function, \ldots may contain extra parameters,
#'   if \code{lambda} is empty, these are parameters passed to the \bold{spatstat}-function \code{\link{density.ppp}}.
#' @param normpower an integer between 0 and 2. If \code{normpower} > 0, the intensity is normalized, see the Details.
#' @return A point pattern (object of class \code{\link{"ppp"}}), with 
#'   additional attribute \code{"htypes"} set to contain the character \code{"w"}.
#' @export
#' @details
#'   If  \code{lambda} is missing, the function will check if the pattern \code{X} 
#'   has previously been marked with an inverse scale factor (by function \code{\link{rescale}}). In that case,
#'   the intensity is calculated as the square inverse scale factor.
#'   
#'   If  \code{lambda} is empty, and \code{X} has no inverse scale factor marks, the intensity
#'   is estimated using the \bold{spatstat}-function \code{\link{density.ppp}}, and extra arguments \ldots
#'   are passed to \code{density.ppp}.
#'    
#'   If  \code{lambda} is given as a function, extra parameters may be passed as \ldots.
#'   
#'   If \code{lambda} is given as a single number, the point pattern gets marked with 
#'   constant intensity, which effectively means that it will be analysed as a homogeneous 
#'   point process with known intensity.
#'   
#'   If  \code{normpower} > 0, the intensity is renormalized, so that \code{\link{Khidden}} yields similar results as
#'   the \bold{spatstat}-function \code{\link{Kinhom}}. The intensity values are then multiplied by 
#'   \deqn{c^{normpower/2}} {c^(normpower/2)}, where 
#'   \deqn{c = area(W)/sum_i(1/\lambda(x_i))}{c = area(W)/sum[i](1/lambda(x[i]))}.
#' @seealso \code{\link{rescaled.ppp}}, \code{\link{retransformed}}    
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

reweighted <- function (X, lambda = NULL, ..., normpower = 0)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if(!is.null(X$marks)) marx <- as.data.frame(X$marks) else marx <- NULL
  if (npts > 0){
    if (is.null(lambda) && !is.null(marx$invscale)) lambda <- marx$invscale^2
    la <- getIntensity(X, lambda, ...,  normpower = normpower)
    if (!is.null (marx$lambda)) marx$lambda <- la 
    else { 
      X$marks <-  as.data.frame(cbind(marx, lambda = la))
      X$markformat <- class(X$marks)
    }
  }
  htypes <- attr(X, "htypes")
  if (!("w" %in% htypes)) attr(X, "htypes") <- c(htypes, "w") 
  return(X)
}

#' Make point pattern locally rescaled 2nd-order stationary
#' 
#' Add mark information for locally rescaled second-order stationary point processes.
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
#' @param normpower an integer between 0 and 2. If \code{normpower} > 0 and \code{invscale} is missing, 
#'   the intensity \code{lambda} is normalized,
#'   see the Details.
#' @return A point pattern (object of class \code{\link{"ppp"}}), with 
#'   additional attribute \code{"htypes"} set to contain the character \code{"s"}.
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
#' @seealso \code{\link{reweighted}}, \code{\link{retransformed}}    
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}


rescaled <- function (X,  invscale = NULL, lambda = NULL, ..., normpower = 0)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if(!is.null(X$marks)) marx <- as.data.frame(X$marks) else marx <- NULL
  if (npts > 0){
  if (is.null(invscale))  {
    if (is.null(lambda)){
       if (!is.null(marx$lambda)) la <- marx$lambda
       else la <- getIntensity(X, NULL, ..., normpower = normpower)
    }  
    else la <- getIntensity(X, lambda, ..., normpower = normpower)   
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
  else{
    X$marks <-  as.data.frame(cbind(marx, invscale = iscale))
    X$markformat <- class(X$marks)
  }
  }
  htypes <- attr(X, "htypes")
  if (!("s" %in% htypes)) attr(X, "htypes") <- c(htypes, "s") 
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
#' @param \ldots optional extra parameters for \code{trafo}, if this is given as a \code{function}.
#' @return A point pattern (object of class \code{\link{"ppp"}}), with 
#'   additional attribute \code{"htypes"} set to contain the character \code{"t"}.
#' @seealso \code{\link{rescaled}}, \code{\link{reweighted}}
#' @details
#'   It is assumed that the transformation is one-to-one, in the sense that the observation window is mapped to itself.
#'   Whether this is true is \emph{not} checked. If you need to analyse transformations that are not one-to-one,
#'   this can be accomplished manually by applying the function \code{\link{coordTransform}} to the pattern \code{X},
#'   and analysing the resulting pattern as a homogeneous point pattern.
#'   
#'   If \code{trafo = "gradx"} or \code{trafo = "grady"}, the backtransformation is obtained from the 
#'   data by inverse cdf transform, under
#'   the assumption that the transformation function affects only the \eqn{x}- or \eqn{y}-coordinate, respectively.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

retransformed <- function (X, trafo = identxy, ...)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if (npts > 0){
  if (is.function(trafo))
  {
    xy <- trafo(X$x, X$y)
    x0 <- xy$x
    y0 <- xy$y
  } 
  else if (trafo %in% c("gradx", "grady"))
    { # just an ordinary inverse cdf transformation, could be more sophisticated at some point
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
  }
  else stop(paste(sQuote("trafo"),
                  "should be a function, or one of the character strings",
                  dQuote("gradx"),"or", dQuote("grady")))
  
  if(!is.null(X$marks)) old.marx <- as.data.frame(X$marks) else old.marx <- NULL
  # gives a funny name for old marks ;-)
  if (!is.null (old.marx$x0)) old.marx$x0 <- x0
  else old.marx <-  as.data.frame(cbind(old.marx, x0 = x0))
  if (!is.null (old.marx$y0)) old.marx$y0 <- y0
  else old.marx <-  as.data.frame(cbind(old.marx, y0 = y0))
  X$marks <- old.marx
  X$markformat <- class(old.marx)
  }
  htypes <- attr(X, "htypes")
  if (!("t" %in% htypes)) attr(X, "htypes") <- c(htypes, "t") 
  return(X)
}  


#' Extract backtransformed pattern
#' 
#' Extract the backtransformed template from a retransformed point process.
#' 
#' @param X the (inhomogeneous) point pattern. Needs to have marks \code{x0} and \code{y0}.
#' @return a point pattern of \bold{spatstat}-class \code{\link{ppp}}
#'  with points in \code{x0} and \code{y0} and marks stripped by these two columns.
#'  The result has an additional attribute \code{"htypes"} set to contain the character \code{"h"},
#'  and is marked with constant intensity and constant inverse scale factor.
#' @export
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' bronzetra <- retransformed(bronzefilter, "gradx")
#' bronzetemplate <- extractRetransformed(bronzetra)
#' plot(bronzetemplate)
extractRetransformed <- function(X)
{
  marx <- marks(X)
  stopifnot (!is.null (marx$x0), !is.null (marx$y0))
  Y <- ppp(marx$x0, marx$y0, window = X$window)
  old.marx <- marks(X)
  old.marx <- old.marx[, !(names(old.marx) %in% c("x0", "y0", "lambda", "invscale"))]
  npts <- npoints(Y)
  la <- npts / area.owin(Y)
  Y$marks <- as.data.frame(cbind(old.marx, lambda = rep(la, npts), invscale = sqrt(la)))
  # marks(Y) <- old.marx[, !(names(old.marx) %in% c("x0", "y0"))]
  return(Y)
}


# @param X point pattern, of class ppp
# @param lambda optional intensity, number, vector, function or image
# @param ... extra parameters for estimation of intensity 
# @param normpower renormalizing power.
#' @rdname sosspp-internal
#' @keywords internal

getIntensity <- function(X, lambda = NULL, ..., normpower = 0)
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  if(is.null(lambda)) {
    # No intensity data provided
    # Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., at="points", leaveoneout=TRUE)
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
  if (normpower != 0) {
    stopifnot ((1 <= normpower) & (normpower <= 2)) 
    renorm.factor <-  (sum(1 / lambda) / (area.owin(X)))^(normpower / 2) 
    lambda <- lambda * renorm.factor
  }
  return(lambda)
}


# @param X point pattern, of class ppp
# @param type assumed type of second-order stationarity
# @return a list (\code{htype}, \code{marx}) with matched type and mark data frame of \code{X}
#' @rdname sosspp-internal
#' @keywords internal


matchtype <- function (X, type = c("w", "t", "s", "h", "hs"))
{
  verifyclass(X, "ppp")
  marx <- X$marks
  if (!is.data.frame(marx)) stop("Given point pattern is not marked hidden 2nd-order stationary.")
  knowntype <-  any(!is.na(match(type, c("w", "t", "s", "h", "hs"))))
  if (!knowntype) stop ("unknown type od hidden 2nd-order stationarity")
  
  # check if X matches type
  htype <- NULL
  if ("w" %in% type) { 
     if (!is.null(marx$lambda)) htype <- "w" }
  else if ("t" %in% type) {
     if (!is.null(marx$x0) && !is.null(marx$y0)) htype <- "t"}
  else if ("s" %in% type) {
     if (!is.null(marx$invscale) || !is.null(marx$lambda)) htype = "s" 
     if (is.null(marx$invscale)) marx <- as.data.frame(cbind(marx, invscale = sqrt(marx$lambda))) }
  else if ("h" %in% type) {  
     la <- rep(npoints(X) / area.owin(X), npoints(X))
     htype <- "h"
     if (is.null(marx$lambda))  marx <- as.data.frame(cbind(marx, lambda = rep(la, npoints(X)))) }
  else if ("hs" %in% type) {  
      iscale <- sqrt(rep(npoints(X) / area.owin(X), npoints(X)))
      htype <- "hs"
      if (is.null(marx$invscale))  marx <- as.data.frame(cbind(marx, invscale = rep(iscale, npoints(X)))) }   
  else htype <- NULL
  if(is.null(htype)) stop ("given inhomogeneity type does not match point pattern")
  return(list(htype = htype, marx = marx))
}



# @param X point pattern, of class ppp, with attribute htype
# @return a list (\code{htype}, \code{marx}) with matched type and mark data frame of \code{X}
# if several types are present, the last one is picked
#' @rdname sosspp-internal
#' @keywords internal

getlasttype <- function (X)
{
  verifyclass(X, "ppp")
  htype <- attr(X, "htypes")
  if (is.vector(htype)) htype <- htype[length(htype)] 
  return(matchtype(X, htype))
}