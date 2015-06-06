# basic methods for sostyppp, like print, extract and replace

#' Check whether an object is a second-order stationarity typed point pattern
#'
#' Checks if an object belongs to class \code{"sostyppp"}.
#'
#' @param x any \code{R} object
#' @return \code{TRUE} if \code{x} belongs to class \code{"sostyppp"}, otherwise \code{FALSE}.
#' @export
#' @seealso \code{\link{sostyppp}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

is.sostyppp <- function(x) inherits(x, "sostyppp")


#' Extract or replace subset of a second-order stationarity typed point pattern
#'
#' An analogon to extraction an replacement of points in an ordinary
#' point pattern, i.e. a spatstat-object of class \code{ppp}.
#'
#' @rdname Extract.sostyppp
# @S3method [ sostyppp
#' @method [ sostyppp
#' @export
#' @param x a sos-typed point pattern, object of class \code{"sostyppp"}.
#' @param i subset index, see spatstat \code{\link[spatstat]{[<-.ppp}}.
# @param j,drop ignored.
#' @seealso \code{\link{sostyppp}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}


"[.sostyppp" <- function(x, i) #, j, drop, ...)
  {
    sostinfo <- attr(x, "sostinfo")
    need.to.rescue.tmarks <- (!is.null(sostinfo$tmarks))
    # attach typemarks to marks
    if (need.to.rescue.tmarks)
    {
      marx <- marks(x)
      tmarx <- as.matrix(sostinfo$tmarks)
      tmarknames <- names(sostinfo$tmarks)

      if (is.null(marx)) {
        marx <- tmarx
        male <- 0
      } else {
        marx <- as.matrix(marx)
        male <- dim(marx)[2]
        marx <- cbind(marx, tmarx)
      }
      x$marks <- as.data.frame(marx)
      y <- NextMethod()
      y <- as.sostyppp(y, "none")
      # now split up the marks again
      marx <- as.matrix(y$marks)
      if (male > 0) {
        # new side effects in spatstat make this a dangerous thing: y looses its type!
        # marks(y) <- marx[, (1:male)] # has to come first, function marks will destroy typemarks!
        y$marks <- marx[, (1:male)]
        sostinfo$tmarks <- as.data.frame(marx[, -(1 : male)], row.names = NULL)
      }
      else {
        y$marks <- NULL
        y$markformat <- "none"
        sostinfo$tmarks <- as.data.frame(marx, row.names = NULL)
      }
      names(sostinfo$tmarks) <- tmarknames 
    }
    else
    {
      y <- NextMethod()
      y <- as.sostyppp(y, "none")
    }
    attr(y, "sostype") <- attr(x, "sostype")
    attr(y, "sostinfo") <- sostinfo
    attr(y, "extra") <- attr(x, "extra")
    return(y)
  }

#' @rdname Extract.sostyppp
# @S3method [<- sostyppp
#'@usage \method{[}{sostyppp} (x, i) <- value
#' @method [<- sostyppp
#' @export
#' @param value Replacement for the subset, a sos-typed point pattern of same type.
#' @details Replacing a subset in a gradient retransformed point patterns does not make much sense
#' since the transformation was calculated from the original data.
"[<-.sostyppp" <-
  function(x, i, value) #j, value)
{
    stopifnot(is.sostyppp(value))

    if(missing(i)) # && missing(j))
      return(value)

    #check if both have same sos-type
    stopifnot(hasType(x, currentType(value)))

    sostyp <- attr(x, "sostype")
    xtras <- attr(x, "extra")
    sostinfo <- attr(x, "sostinfo")
    tmarknames <- names(sostinfo$tmarks)

    marx <- marks(x)
    marknames <- names(marx)

    need.to.rescue.tmarks <- (!is.null(sostinfo$tmarks))
    if (need.to.rescue.tmarks)
    {
      #attach typemarks to marks in order to use spatstats replace mechanism
      tmarx <- as.matrix(sostinfo$tmarks)
       if (is.null(marx)) {
          marx <- tmarx
          male <- 0
      } else {
        marx <- as.matrix(marx)
        male <- dim(marx)[2]
        marx <- cbind(marx, tmarx)
      }

      names(marx) <- NULL
      mary <- as.matrix(marks(value))
      maly <- dim(mary)[2]
      tvmarx <- as.matrix(attr(value, "sostinfo")$tmarks)
      if (maly>0) mary <- cbind(mary, tvmarx)  else mary <- tvmarx
      names(mary) <- NULL

      # marks function is not inheritable - securing valuables before using it...

      marks(x) <- as.data.frame(marx, row.names = NULL)
      marks(value) <- as.data.frame(mary, row.names = NULL)
      y <- NextMethod()
      y <- as.sostyppp(y, "none")
      
      marx <- as.matrix(y$marks)
      if (male > 0) {
        y$marks <- as.data.frame(marx[, (1:male)], row.names = NULL) # has to come first, function marks will destroy typemarks!
        sostinfo$tmarks <- as.data.frame(marx[, -(1 : male)], row.names = NULL)
      }
      else {
        y$marks <- NULL
        sostinfo$tmarks <- as.data.frame(marx, row.names = NULL)
      }
      names(sostinfo$tmarks) <- tmarknames
    }
    else
    {
      y <- NextMethod()
      y <- as.sostyppp(y)
    }
    if(!is.null(y$marks)) names(y$marks) <- marknames
    attr(y, "sostinfo") <- sostinfo
    attr(y, "sostype") <- sostyp
    attr(y, "extra") <- xtras
    invisible(y)
}

#' Print brief details of a second-order stationarity typed point pattern
#'
#' Extends spatstat method \code{\link[spatstat]{print.ppp}} by type of second-order stationarity.
#'
#' @param x sos-typed point pattern.
#' @param ... ignored
# @S3method print sostyppp
#' @method print sostyppp
#' @export
#' @seealso \code{\link{print.ppp}} for the print method of class ancestor \code{ppp}.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

print.sostyppp <- function(x, ...)
{
  print.ppp(x)
  sost <- attr(x, "sostype")
  if(length (sost) > 0)
    cat("pattern is", longTypeName(sost[1]),"second-order stationary","\n")
  else cat("pattern is second-order stationary of unassigned type","\n")
  further <- sost[-1]
  if (length(further > 0)) 
    cat("additional types:",sapply(further, longTypeName),"\n")
}

#' Manipulate s.o.s. type information
#'
#' Access the encrypted type information of second-order stationarity typed point
#' patterns, objects of class \code{sostyppp}.
#' @param x point pattern, of class sostyppp
#' @param type character, type of second-order stationarity. See \code{\link{sostyppp}}
#' for details.
# @return logical
#' @export
# @rdname sostyppp-types
# @alias sos-type functions
# @keywords internal
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

hasType <- function (x, type = .TYPES)
{
  knowntype <-  any(!is.na(match(type, .TYPES)))
  if (!knowntype) stop ("unknown type of hidden 2nd-order stationarity")
  #return(type %in% .gettype(attr(x, "sostype"))$all)
  return(type %in% (attr(x, "sostype")))
}


# @param x point pattern, of class sostyppp
# @return character, giving the last second-order stationarity type X was assigned to
# if several types are present, the last one is picked
# @rdname sostyppp-types
# @keywords internal
#' @export
# @alias sos-type functions
#' @rdname hasType

currentType <- function (x)
{
  #return(.gettype(attr(x, "sostype"))$last)
  return(attr(x, "sostype")[1])
}

#@param value type character
#@usage currentType(x) <- value
#'@export
#'@rdname sostatpp-internal
#'@usage currentType(x) <- value
#'@keywords internal

"currentType<-" <- function(x, value)
{
  knowntype <-  any(!is.na(match(value, .TYPES)))
  if (!knowntype) stop ("unknown type of hidden 2nd-order stationarity")
  types <- attr(x, "sostype")
  if("none" %in% types) types <- types[-match("none", types)]
  if(value %in% types) types <- types[-match(value, types)]
  attr(x, "sostype") <- c(value, types)
  x
}  

# @param typecode integer. encrypted information on last type, and all types
# @return a list (\code{last}, \code{all}) of character vectors giving the type
# if several types are present, the last one is picked
# @rdname sostatpp-internal
# @keywords internal

.TYPENAMES  <- c("reweighted", "retransformed", "rescaled", "homogeneous", 
  "scaled-homogeneous", "not specified")
.TYPES  <- c("w", "t", "s", "h", "hs", "none")

# return long name
#
#' @rdname sostatpp-internal
#' @keywords internal

longTypeName <- function(type){
  typeno <- match(type, .TYPES)
  if (!is.na(typeno)) .TYPENAMES[typeno] else ""
}

# @param typecode character
# @param type character, giving the type to be entered in encrypted info
# @return integer, updated encrypted type information
# if several types are present, the last one is picked
# @export
#' @rdname sostatpp-internal
#' @keywords internal

.settype <- function(type, typecode)
{
  knowntype <-  any(!is.na(match(type, .TYPES)))
  if (!knowntype) stop ("unknown type of hidden 2nd-order stationarity")
  if("none" %in% typecode) typecode <- typecode[-match("none", typecode)]
  if(type %in% typecode) typecode <- typecode[-match(type, typecode)]
  c(type, typecode)
}  

# .settype_oldsystem <- function(type, typecode)
# {
#   tindex <- which(.TYPES == type)
#   stopifnot(length(tindex) == 1)
#   if (length(typecode) < 1) { # new typecode
#     return(.TYBITS[tindex] * (.TYLAST+1))
#   }
#   else {
#     all <- typecode %% .TYLAST
#     contained <- (all %/% .TYBITS) %% 2 == 1
#     if (!contained[tindex]) all <- all + .TYBITS[tindex]
#     return (.TYBITS[tindex] * .TYLAST + all)
#   }
# }



# documentation of class sostyppp -------------------------------------------

#'@name sostyppp
#'@aliases sostyppp
#'@title Class for Second-Order Stationarity-Typed Planar Point Patterns
#'@description
#'  A class \code{"sostyppp"} to represent a two-dimensional point
#'  pattern, with extra information about hidden second-order
#'  stationarity. Built on top of the spatstat-class \code{"\link{ppp}"}.
#'@details
#'  This class represents a two-dimensional point pattern dataset with
#'  extra information relevant for analysis as a realization of a
#'  second-order stationary point process.
#'
#'  From its ancestor, spatstat-class \code{"\link[spatstat]{ppp}"}, an object of type
#'  \code{sostyppp} inherits the elements
#'  \tabular{ll}{
#'    \code{x} \tab vector of \eqn{x} coordinates of data points
#'    \cr\code{y} \tab vector of \eqn{y} coordinates of data points
#'    \cr\code{n} \tab number of points
#'    \cr\code{window} \tab window of observation (an object of class \code{\link{owin}})
#'    \cr\code{marks} \tab vector or data frame of marks
#'  }
#'  Additionally, it has attributes
#' \tabular{ll}{
#'    \code{sostype} \tab integer, encrypts the type of stationarity that is
#'    assumed in analysis, see below
#'    \cr\code{sostinfo} \tab a list containing information relevant to
#'    for calculating second-order statistics,
#'    \cr\tab such as the intensity function at  the points of the pattern,
#'    \cr\code{extra} \tab a container (list) of additional information, intended for internal use.
#'    \cr\tab This container is copied in subsetting operations.\cr
#'  }
#'
#'  The type of second-order stationarity of an object \code{X} of class \code{sostyppp}
#'  is returned as \code{character} by the function \code{\link{currentType}}. Whether \code{X} has a given type of
#'  second-order stationarity, can be checked with \code{\link{hasType}}.
#'
#'  Possible types of second-order stationarity (s.o.s.) are
#'  \tabular{ll}{
#'    \code{"w"}  \tab intensity reweighted s.o.s.
#'                \cr\tab \code{tinfo} contains a data frame \code{tmarks},
#'                 with a variable \code{intens} giving the 
#'               \cr\tab intensity evaluated in each data point.
#'    \cr\code{"t"} \tab obtained by coordinate transformation
#'               \cr\tab attribute \code{sostinfo} contains a function 
#'               \code{backtrafo(x,y)} that yields the backtransformed 
#'               \cr\tab pattern, and a character variable \code{gradient}, taking 
#'               values \code{"gradx"} or \code{"grady"}
#'                \cr\tab if the transformation depends only on one coordinate.
#'    \cr\code{"s"} \tab locally rescaled s.o.s.
#'                  \cr\tab \code{sostinfo} contains a data frame \code{tmarks}, with a variable
#'                    \code{invscale}, 
#'                    \cr\tab the inverse scale factor in each data point\cr
#'    \cr\code{"h"}  \tab homogeneous, to be evaluated with standard methods
#'    \cr\code{"hs"} \tab homogeneous, to be evaluated with scale invariant statistics.
#'    }
#'  An object of class \code{sostyppp} can carry type information for several types simultaneously.
#'  Thus, function \code{\link{hasType}} may return \code{TRUE} for different arguments.
#'
#'@author Ute Hahn  \email{ute@@imf.au.dk}
#'
#'@keywords spatial attribute
#'
NA
