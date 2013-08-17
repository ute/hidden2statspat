# Funktions for coordinate transformation of point processes and windows

# transformation has to be given as a function on a data frame or list with 
# elements x and y

# helper functions, transform image, subdivide polygons----------------------------------

#' internal functions of hidden2ndorder package
#' @name internal_functions
#' @description
#' functions for use by the package's functions, partly slightly 
#' documented in the source code
# @param xy list or dataframe with elements \code{x} and \code{y}
#' @rdname hidden2ndorder-internal
#' @keywords internal

identxy <- function (xy) return (xy)

#' @rdname hidden2ndorder-internal
#' @keywords internal
mapstructxy <- function(X, mapxy = identxy)
{
  newxy <- mapxy(X)
  X$x <- newxy$x
  X$y <- newxy$y
  return(X)
}

# @param poly list or dataframe with elements \code{x} and \code{y}
#' @rdname hidden2ndorder-internal
#' @keywords internal

peripolyxy <- function (poly)
{
  xx <- c(poly$x, poly$x[1])
  yy <- c(poly$y, poly$y[1])
  return(sum(sqrt(diff(xx)^2 + diff(yy)^2)))
}

# @param poly list or dataframe with elements \code{x} and \code{y}
# @param newlen maximal length of edges in refined polygon
#' @rdname hidden2ndorder-internal
#' @keywords internal

subdivpolyxy <- function (poly, newlen)
{
  addi <- function(x, dx, n) x + (0 : (n-1)) * dx / n
  xx <- poly$x
  yy <- poly$y
  dx <- diff(c(xx, xx[1]))
  dy <- diff(c(yy, yy[1]))
  lengths <- sqrt(dx^2 + dy^2)
  ndiv <- ceiling(lengths / newlen) 
  xnew <- unlist(mapply(addi, xx, dx, ndiv))
  ynew <- unlist(mapply(addi, yy, dy, ndiv))
  poly$x <- xnew
  poly$y <- ynew
  poly$perimeter <- sum(lengths)
  return(poly)
}

#' Refine a Polygonal Window
#'
#' Generate a finer polygonal window by subdividing edges.
#'
#' @param X A polygonal window (object of class \code{"owin"} and of type \code{"polygonal"}).
#' @param edgelen The maximum edge length (see `Details') of the resulting 
#' polygon, defaults to perimeter/50  of the largest subpolygon.
#' @details Every edge is subdivided into segments of same length, with length 
#' less or equal to \code{edgelen}.
#' @export
#' @keywords manip
#' @keywords spatial
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' # refine the polygonal window letterR from the spatstat package
#' fineR <- refinepoly(letterR, edgelen = 0.1)
#' summary(letterR)
#' plot(letterR)
#' summary(fineR)
#' plot(fineR)

refinepoly <- function (X, edgelen = NULL) 
{
  verifyclass(X, "owin")
  poly <- as.polygonal(X)
  if (is.null(edgelen)) edgelen <- max(sapply(poly$bdry, peripolyxy)) / 50
  newbdry <- lapply(poly$bdry, subdivpolyxy, newlen = edgelen)
  poly$bdry <- newbdry
  return(poly)
}

#' Apply Coordinate Transformation
#' 
#' Apply a coordinate transformation to a point pattern or a window  
#' 
#'
#' @param X Object to be transformed, currently of class \code{"\link{im}"}, 
#'        \code{"\link{owin}"} or \code{"\link{ppp}"}.
#' @param \ldots arguments passed to class methods.
#' @export
#' @details The functions \code{trafoxy} (and \code{invtrafoxy}) take one argument
#'      that is a list or a data frame with elements \code{x} and \code{y} and return
#'      an argument of same type.
#' @seealso \code{\link{coordTransform.im}}, \code{\link{coordTransform.owin}},
#' \code{\link{coordTransform.ppp}}

coordTransform <- function(X,  ...)
{
  UseMethod("coordTransform", X,  ...)
}


#' Apply Coordinate Transformation to a Pixel Image
#'
#' Apply any coordinate transformation to
#' a pixel image. Generalizes the spatstat function 
#' \code{\link{affine.im}}. 
#'
#' @param X Pixel image (object of class \code{"im"}).
#' @param trafoxy the coordinate transformation function, defaults to the identical 
#'  map \code{\link{identxy}}, see `Details'.
#' @param invtrafoxy the inverse to \code{trafoxy}, defaults to the identical 
#'  map \code{\link{identxy}}.
#' @param ... Optional arguments passed to \code{\link{as.mask}} controlling the 
#'  pixel resolution of the transformed image. 
#' @details The functions \code{trafoxy} (and \code{invtrafoxy}) take one argument
#'      that is a list or a data frame with elements \code{x} and \code{y} and return
#'      an argument of same type.
#' @keywords spatial
#' @keywords manip
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @S3method coordTransform im
#' @examples
#' # transformation of the polygonal window letterR from the spatstat package
#' # as image
#' trafo <- function(xy) list(x = xy$x^1.5 + 0.5*xy$y, y = xy$y*.5)
#' invtrafo <- function(xy) list(x = (xy$x - xy$y) ^(1/1.5), y = xy$y*2)
#' imR <- as.im(letterR)
#' plot(imR)
#' plot(coordTransform.im(imR, trafo, invtrafo))

coordTransform.im <- function (X, trafoxy = identxy, invtrafoxy = identxy, ...) 
{
  verifyclass(X, "im")
  {
    newframe <- bounding.box.xy(trafoxy(corners(X)))
    W <- if (length(list(...)) > 0) 
      as.mask(newframe, ...)
    else as.mask(newframe, eps = with(X, min(xstep, ystep)))
    naval <- switch(X$type, factor = {factor(NA, levels = levels(X))}, 
      integer = NA_integer_, logical = as.logical(NA_integer_), 
      real = NA_real_, complex = NA_complex_, character = NA_character_, 
      NA)
    Y <- as.im(W, value = naval)
    xx <- as.vector(rasterx.im(Y))
    yy <- as.vector(rastery.im(Y))
    pre <- invtrafoxy(list(x = xx, y = yy))
    Y$v[] <- lookup.im(X, pre$x, pre$y, naok = TRUE)
    return(Y)
  }
}  


# transform owin ----------------------------------------------------------

#' Apply Coordinate Transformation to Window
#'
#' Apply any coordinate transformation to a window.
#' Compare the spatstat function \code{\link{affine.owin}}. 
#'
#' @param X Window (object of class \code{"owin"}).
#' @param trafoxy The coordinate transformation function, defaults to the identical 
#'        map \code{\link{identxy}}, see also `Details'.
#' @param invtrafoxy The inverse to \code{trafoxy}, defaults to the identical map 
#'       \code{\link{identxy}}. Only needed if \code{X} is a pixel image.
#' @param isAffine If true, no subdivision of rectangles or polygons, see the `Details'.
#' @param \ldots Optional arguments passed to \code{\link{as.mask}} controlling the 
#'      pixel resolution of the transformed window, if \code{X} is a binary pixel 
#'      mask, or to \code{\link{subdivpoly}}, if \code{X} is a polygon or a 
#'      rectangle, see the `Details'. 
#' @details The functions \code{trafoxy} (and \code{invtrafoxy}) take one argument
#'      that is a list or a data frame with elements \code{x} and \code{y} and return
#'      an argument of same type.
#'      
#'      If the window is a rectangle or polygon, it is converted into a 
#'      polygon which is subsequently refined using \code{\link{subdivpoly}}, 
#'      unless \code{isAffine=TRUE}. This is only necessary if the coordinate 
#'      transform is not affine, in order to achieve a better approximation of 
#'      the transformed window. Note however that \code{coordTransform.owin} does 
#'      not check whether \code{trafoxy} really is affine. The option 
#'      \code{isAffine=TRUE} may therefore also be used to save memory and 
#'      computation time when transforming xy-rectangles with a transformation 
#'      that preserves axe-parallel rectangles.
#' @seealso \code{\link{coordTransform.im}}, which is called if \code{X} is a pixel image. 
#' @export
#' @S3method coordTransform owin
#' @keywords manip
#' @keywords spatial
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' # transformation of the polygonal window letterR from spatstat
#' 
#' # a nonlinear transformation 
#' trafo <- function(xy) list(x = xy$x^1.5 + 0.5*xy$y, y = xy$y*.5)
#' 
#' # transform letterR and mark the boundary corners, 
#' # dummy variable for silent return
#' plot(letterR)
#' dummy <- lapply(letterR$bdry, points, col = "red")
#' mappedR <- coordTransform.owin(letterR, trafo, edgelen = 0.1)
#' plot(mappedR)
#' dummy <- lapply(mappedR$bdry, points, col = "red")
#' 
#' # no refinement, assuming that trafo is affine
#
#' mappedRcoarse <- coordTransform.owin(letterR, trafo, isAffine = T)
#' plot(mappedRcoarse, add=T)
#' dummy <- lapply(mappedRcoarse$bdry, points, col = "green")
#' # not much of a difference, though...
#' 
#' # now as mask image. We need the inverse transformation, too
#' 
#' invtrafo <- function(xy) list(x = (xy$x - xy$y) ^(1/1.5), y = xy$y*2)
#' wimR <- as.owin(as.im(letterR))
#' mappedimR <- coordTransform.owin(wimR, trafo, invtrafo)
#' plot(mappedimR) 
#' # compare with the polygonal window from before
#' plot(mappedR, col = "green", add = T)


coordTransform.owin <- function (X, trafoxy = identxy, invtrafoxy = NULL, isAffine = F, ...)
{
    verifyclass(X, "owin")
    if (X$type %in% c("rectangle", "polygonal"))
      {
        P <- if(isAffine) as.polygonal(X) else refinepoly(X, ...)
        newbdry <- lapply(P$bdry, mapstructxy, mapxy = trafoxy)
        # recalculate areas
        newerbdry <- lapply(newbdry, 
                function(b) {b$area <-NULL; b$perimeter <- NULL;
                         b$area <- area.xypolygon(b); return(b)})  
        P$bdry <- newerbdry
        P$xrange <- range(sapply(P$bdry, function(p) p$x))
        P$yrange <- range(sapply(P$bdry, function(p) p$y))
        if(isAffine & (X$type == "rectangle")) P <- as.rectangle(P)
        return(P)
       }  else if (X$type == "mask") {
         stopifnot(!is.null(invtrafoxy))
         return(as.owin(coordTransform.im(as.im(X), trafoxy, invtrafoxy, ...))) 
       } else stop("Unrecognised window type")
}


# transform ppp ----------------------------------------------------------

#' Apply Coordinate Transformation to Point Pattern
#'
#' Apply any coordinate transformation to a point pattern.
#' Compare the spatstat function \code{\link{affine.ppp}}. 
#'
#' @param X a point pattern (object of class \code{"ppp"}).
#' @param trafoxy the coordinate transformation function, defaults to the identical 
#'        map \code{\link{identxy}}, see also `Details'.
#' @param \ldots Optional arguments passed to \code{\link{coordTransform.owin}}. 
#' @details The functions \code{trafoxy} takes one argument, namely
#'      a list or a data frame with elements \code{x} and \code{y}, and returns
#'      an argument of same type.
#'      
#'      If the window of $X$ is a pixel mask, an additional inverse transformation
#'      has to be provided, see \code{\link{coordTransform.owin}}.
#' @seealso \code{\link{coordTransform.owin}} for the transformation of windows     
#' @export
#' @S3method coordTransform ppp
#' @keywords manip
#' @keywords spatial
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @examples
#' # transforming a Poisson point pattern by inhomogeneous scaling 
#' # around origo
#' 
#' galaxytrafo <- function(xy) { 
#'   xx <- xy$x; yy <- xy$y; dd <- xx^2+yy^2
#'   return(list(x = dd*xx, y = dd*yy))
#'   }
#'
#' pp <- rpoispp(100, win = shift(square(2), c(-1,-1)))
#' ppstar <- coordTransform.ppp(pp, galaxytrafo)
#' plot(pp, pch=16, cex=.5)
#' plot(ppstar, pch=16, cex=.5)


coordTransform.ppp <- function(X, trafoxy = identxy,  ...)#invtrafoxy=NULL, isAffine = F, ...)
{
  Pnew <- mapstructxy(X, trafoxy)
  Pnew$window <- coordTransform.owin(X$window, trafoxy,  ...)#invtrafoxy, isAffine, ...)
  return(Pnew)
}