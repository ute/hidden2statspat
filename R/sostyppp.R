# basic methods for sostyppp, like print, extract and replace

#' Check whether an object is a second-order stationarity typed point pattern
#'
#' Checks if an object belongs to class \code{"sostyppp"}.
#'
#' @param x any \code{R} object
#' @return \code{TRUE} if \code{x} belongs to class \code{"sostyppp"}, otherwise \code{FALSE}.
#' @export
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

is.sostyppp <- function(x) inherits(x, "sostyppp")


#' Extract or replace subset of a second-order stationarity typed point pattern
#'
#' An analogon to extraction an replacement of points in an ordinary
#' point pattern, i.e. a spatstat-object of class \code{ppp}.
#'
#' @rdname Extract.sostyppp
#' @S3method [ sostyppp
#' @method [ sostyppp
# @export
#' @param x a sos-typed point pattern, object of class \code{"sostyppp"}.
#' @param i subset index, see spatstat \code{\link[spatstat]{[<-.ppp}}.
# @param j,drop ignored.
#' @seealso \code{\link{sostyppp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}


"[.sostyppp" <- function(x, i) #, j, drop, ...)
  {
    need.to.rescue.tmarks <- (!is.null(x$sostinfo$tmarks))
    # attach typemarks to marks
    if (need.to.rescue.tmarks)
    {
      marx <- marks(x)
      tmarx <- as.matrix(x$sostinfo$tmarks)
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
      y$sostinfo <- x$sostinfo
      # now split up the marks again
      marx <- as.matrix(y$marks)
      if (male > 0) {
        # new side effects in spatstat make this a dangerous thing: y looses its type!
        # marks(y) <- marx[, (1:male)] # has to come first, function marks will destroy typemarks!
        y$marks <- marx[, (1:male)]
        y$sostinfo$tmarks <- as.data.frame(marx[, -(1 : male)], row.names = NULL)
      }
      else {
        y$marks <- NULL
        y$markformat <- "none"
        y$sostinfo$tmarks <- as.data.frame(marx, row.names = NULL)
      }
      names(y$sostinfo$tmarks) <- names(x$sostinfo$tmarks)
    }
    else
    {
      y <- NextMethod()
      y <- as.sostyppp(y, "none")
      y$sostinfo <- x$sostinfo
    }
    y$sostype <- x$sostype
    y$extra <- x$extra
    return(y)
  }

#' @rdname Extract.sostyppp
#' @S3method [<- sostyppp
#'@usage \method{[}{sostyppp} (x, i) <- value
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
    stopifnot(has.type(x, currenttype(value)))

    sostyp <- x$sostype
    xtras <- x$extra
    sosinfo <- x$sostinfo
    tmarknames <- names(x$sostinfo$tmarks)

    marx <- marks(x)
    marknames <- names(marx)

    need.to.rescue.tmarks <- (!is.null(x$sostinfo$tmarks))
    if (need.to.rescue.tmarks)
    {
      #attach typemarks to marks in order to use spatstats replace mechanism
      tmarx <- as.matrix(x$sostinfo$tmarks)
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
      tvmarx <- as.matrix(value$sostinfo$tmarks)
      if (maly>0) mary <- cbind(mary, tvmarx)  else mary <- tvmarx
      names(mary) <- NULL

      # marks function is not inheritable - securing valuables before using it...

      marks(x) <- as.data.frame(marx, row.names = NULL)
      marks(value) <- as.data.frame(mary, row.names = NULL)
      y <- NextMethod()
      y <- as.sostyppp(y, "none")
      y$sostinfo <- sosinfo

      marx <- as.matrix(y$marks)
      if (male > 0) {
        y$marks <- as.data.frame(marx[, (1:male)], row.names = NULL) # has to come first, function marks will destroy typemarks!
        y$sostinfo$tmarks <- as.data.frame(marx[, -(1 : male)], row.names = NULL)
      }
      else {
        y$marks <- NULL
        y$sostinfo$tmarks <- as.data.frame(marx, row.names = NULL)
      }
    }
    else
    {
      y <- NextMethod()
      y <- as.sostyppp(y)
      y$sostinfo <- x$sostinfo
    }
    if(!is.null(y$marks)) names(y$marks) <- marknames
    y$sostype <- sostyp
    y$extra <- xtras
    names(y$sostinfo$tmarks) <- tmarknames
    invisible(y)
}

#' Print brief details of a second-order stationarity typed point pattern
#'
#' Extends spatstat method \code{\link[spatstat]{print.ppp}} by type of second-order stationarity.
#'
#' @param x sos-typed point pattern.
#' @param ... ignored
#' @S3method print sostyppp
#' @method print sostyppp
# @export
#' @seealso \code{\link{print.ppp}} for the print method of class ancestor \code{ppp}.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

print.sostyppp <- function(x, ...)
{
  print.ppp(x)
  if(length (currenttypeno(x)) > 0)
    cat("pattern is",.TYPENAMES[currenttypeno(x)],"second-order stationary","\n")
  else cat("pattern is second-order stationary of unassigned type","\n")
  further <- furthertypeno(x)
  if (length(further>0)) cat("additional types:",.TYPENAMES[further],"\n")
}

#' Manipulate s.o.s. type information
#'
#' Access the encrypted type information of second-order stationarity typed point
#' patterns, objects of class \code{sostyppp}.
#' @param x point pattern, of class sostyppp
#' @param type character, type of second-order stationarity. See \code{\link{sostyppp.object}}
#' for details.
# @return logical
#' @export
# @rdname sostyppp-types
# @alias sos-type functions
# @keywords internal
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

has.type <- function (x, type = .TYPES)
{
  knowntype <-  any(!is.na(match(type, .TYPES)))
  if (!knowntype) stop ("unknown type of hidden 2nd-order stationarity")
  return(type %in% .gettype(x$sostype)$all)
}


# @param x point pattern, of class sostyppp
# @return character, giving the last second-order stationarity type X was assigned to
# if several types are present, the last one is picked
# @rdname sostyppp-types
# @keywords internal
#' @export
# @alias sos-type functions
#' @rdname has.type

currenttype <- function (x)
{
  return(.gettype(x$sostype)$last)
}

# @return integer, index number of second-order stationarity type
# @rdname sostyppp-types
#' @rdname sostatpp-internal
#' @keywords internal
# @export

currenttypeno <- function (x)
{
  return(.gettype(x$sostype)$lastno)
}


# @param x point pattern, of class sostyppp
# @return vector of type numbers, all but the current type
# if several types are present, the last one is picked
# @rdname sostyppp-types
#' @rdname sostatpp-internal
#' @keywords internal
# @export

furthertypeno <- function (x)
{
  return(.gettype(x$sostype)$furtherno)
}


# hand made bit operations, certainly much too awkward...
# @param typecode integer. encrypted information on last type, and all types
# @return a list (\code{last}, \code{all}) of character vectors giving the type
# if several types are present, the last one is picked
#' @rdname sostatpp-internal
#' @keywords internal
# \code{"ppp"}

.TYPENAMES  <- c("reweighted", "retransformed", "rescaled", "homogeneous", 
  "scaled-homogeneous", "not specified")
.TYPES  <- c("w", "t", "s", "h", "hs", "none")
.TYBITS <- c(1, 2, 4, 8, 16, 32)
.TYLAST <- 64

.gettype <- function (typecode)
{
  last <- typecode %/% .TYLAST
  lindex <- (last %/% .TYBITS) %% 2 == 1
  all <- typecode %% .TYLAST
  contained <- (all %/% .TYBITS) %% 2 == 1
  return(list(last = .TYPES[lindex], all = .TYPES[contained],
              lastno = which(lindex), furtherno = which(!lindex & contained) ))
}

# @param typecode integer, encrypted information on last type, and all types so far
# @param type character, giving the type to be entered in encrypted info
# @return integer, updated encrypted type information
# if several types are present, the last one is picked
# @export
#' @rdname sostatpp-internal
#' @keywords internal

.settype <- function(type, typecode)
{
  tindex <- which(.TYPES == type)
  stopifnot(length(tindex) == 1)
  if (length(typecode) < 1) { # new typecode
    return(.TYBITS[tindex] * (.TYLAST+1))
  }
  else {
    all <- typecode %% .TYLAST
    contained <- (all %/% .TYBITS) %% 2 == 1
    if (!contained[tindex]) all <- all + .TYBITS[tindex]
    return (.TYBITS[tindex] * .TYLAST + all)
  }
}



# documentation of class sostyppp -------------------------------------------

#'@name sostyppp.object
#'@aliases sostyppp.object sostyppp
#'@title Class of Second-Order Stationarity-Typed Point Patterns
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
#'  Additionally, it contains the elements
#' \tabular{ll}{
#'    \code{sostype} \tab integer, encrypts the type of stationarity that is
#'    assumed in analysis, see below
#'    \cr\code{tinfo} \tab a list containing information relevant to
#'    the type of stationarity,
#'    \cr\code{extra} \tab a container (list) of additional information, intended for internal use.
#'    This container is copied in subsetting operations.\cr
#'  }
#'
#'  The type of second-order stationarity of an object \code{X} of class \code{sostyppp}
#'  is returned as \code{character} by the function \code{\link{currenttype}}. Whether \code{X} has a given type of
#'  second-order stationarity, can be checked with \code{\link{has.type}}.
#'
#'  Possible types of second-order stationarity (s.o.s.) are
#'  \tabular{ll}{
#'    \code{"w"}  \tab intensity reweighted s.o.s.
#'                \cr\tab \code{tinfo} contains a data frame \code{tmarks}, with a variable \code{intens} with intensity evaluated
#'               in each data point
#'    \cr\code{"t"} \tab obtained by coordinate transformation
#'               \cr\tab \code{tinfo} contains a function \code{backtrafo(x,y)} that yields the backtransformed pattern.
#'                \cr\tab and a character variable \code{gradient}, taking values \code{"gradx"} or \code{"grady"}
#'                if the transformation depends only on one coordinate.
#'    \cr\code{"s"} \tab locally rescaled s.o.s.
#'                  \cr\tab \code{tinfo} contains a data frame \code{tmarks}, with a variable  \code{invscale}, the inverse scale factor in each data point\cr
#'    \cr\code{"h"}  \tab homogeneous, to be evaluated with standard methods
#'    \cr\code{"hs"} \tab homogeneous, to be evaluated with scale invariant statistics.
#'    }
#'  An object of class \code{sostyppp} can carry type information for several types simultaneously.
#'  Thus, function \code{\link{has.type}} may return \code{TRUE} for different arguments.
#'
#'@author Ute Hahn  \email{ute@@imf.au.dk}
#'
#'@keywords spatial attribute
#'
NA
