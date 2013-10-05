# basic methods for sostpp, like print, extract and replace 

#' Check whether an object is a second-order stationarity typed point pattern
#' 
#' Checks if an object belongs to class \code{"sostpp"}.
#' 
#' @param x any \code{R} object
#' @return \code{TRUE} if \code{x} belongs to class \code{"sostpp"}, otherwise \code{FALSE}.
#' @export
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

is.sostpp <- function(x) inherits(x, "sostpp")


#' Extract or replace subset of a second-order stationarity typed point pattern
#' 
#' An analogon to extraction an replacement of points in an ordinary 
#' point pattern, i.e. a spatstat-object of class \code{ppp}.
#' 
#' @rdname Extract.sostpp
#' @S3method [ sostpp
#' @method [ sostpp
# @export
#' @param x a sos-typed point pattern, object of class \code{"sostpp"}. 
#' @param i subset index, see spatstat \code{\link{[<-.ppp}}.
# @param j,drop ignored.
#' @seealso \code{\link{sostpp.object}} for details on the class.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}


"[.sostpp" <- function(x, i) #, j, drop, ...) 
  {
    # attach typemarks to marks
    marx <- as.data.frame(marks(x))
    male <- dim(marx)[2]
    if (male>0) marx <- cbind(marx, x$typemarks)  else marx <- x$typemarks
    x$marks <- marx
    y <- NextMethod()
    y <- as.sostpp.ppp(y)
    marx <- y$marks
    if (male > 0) {
      marks(y) <- marx[, (1:male)] # has to come first, function marks will destroy typemarks!
      y$typemarks <- as.data.frame(marx[, -(1 : male)]) 
    }
    else {
      marks(y) <- NULL
      y$typemarks <- as.data.frame(marx)
    }
    y$sostype <- x$sostype
    y$extra <- x$extra
    names(y$typemarks) <- names(x$typemarks)
    class(y)<- c("sostpp", class(y))
    return(y)
  }

#' @rdname Extract.sostpp
#' @S3method [<- sostpp
#'@usage \method{[}{sostpp} (x, i) <- value 
#' @export
#' @param value Replacement for the subset, a sos-typed point pattern of same type. 

"[<-.sostpp" <-
  function(x, i, value) #j, value) 
{
    
    stopifnot(is.sostpp(value))
    
    if(missing(i)) # && missing(j))
      return(value)
    
    #check if both have same sos-type
    stopifnot(has.type(x, currenttype(value)))
    
    #attach typemarks to marks
    marx <- as.data.frame(marks(x))
    male <- dim(marx)[2]
    if (male>0) marx <- cbind(marx, x$typemarks)  else marx <- x$typemarks
    mary <- as.data.frame(marks(value))
    maly <- dim(mary)[2]
    if (maly>0) mary <- cbind(mary, value$typemarks)  else mary <- value$typemarks
   
    # marks function is not inheritable - securing valuables before using it...
    typemarknames <- names(x$typemarks)
    sostyp <- x$sostype
    xtras <- x$extra
    
    marks(x) <- marx
    marks(value) <- mary
    y <- NextMethod()
    
    marx <- y$marks
    if (male > 0) {
      marks(y) <- marx[, (1:male)] # has to come first, function marks will destroy typemarks!
      y$typemarks <- as.data.frame(marx[, -(1 : male)]) 
    }
    else {
      marks(y) <- NULL
      y$typemarks <- as.data.frame(marx)    
    }
    y$sostype <- sostyp
    y$extra <- xtras
    names(y$typemarks) <- typemarknames
    class(y)<- c("sostpp", class(y))
    invisible(y)
} 

#' Print brief details of a second-order stationarity typed point pattern
#' 
#' Extends spatstat method \code{\link{print.ppp}} by type of second-order stationarity.
#' 
#' @param x sos-typed point pattern.
#' @param ... ignored
#' @S3method print sostpp
#' @method print sostpp
# @export
#' @seealso \code{\link{print.ppp}} for the print method of class ancestor \code{ppp}.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

print.sostpp <- function(x, ...)
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
#' patterns, objects of class \code{sostpp}.
#' @param x point pattern, of class sostpp
#' @param type character, type of second-order stationarity. See \code{\link{sostpp.object}}
#' for details.
# @return logical
#' @export
# @rdname sostpp-types
# @alias sos-type functions
# @keywords internal
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

has.type <- function (x, type = .TYPES)
{
  knowntype <-  any(!is.na(match(type, .TYPES)))
  if (!knowntype) stop ("unknown type of hidden 2nd-order stationarity")
  return(type %in% .gettype(x$sostype)$all) 
}


# @param x point pattern, of class sostpp
# @return character, giving the last second-order stationarity type X was assigned to
# if several types are present, the last one is picked
# @rdname sostpp-types
# @keywords internal
#' @export
# @alias sos-type functions
#' @rdname has.type

currenttype <- function (x)
{
  return(.gettype(x$sostype)$last) 
}

# @return integer, index number of second-order stationarity type
# @rdname sostpp-types
#' @rdname sostatpp-internal
#' @keywords internal
# @export

currenttypeno <- function (x)
{
  return(.gettype(x$sostype)$lastno) 
}


# @param x point pattern, of class sostpp
# @return vector of type numbers, all but the current type
# if several types are present, the last one is picked
# @rdname sostpp-types
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

.TYPENAMES  <- c("reweighted", "retransformed", "rescaled", "homogeneous", "scaled-homogeneous")
.TYPES  <- c("w", "t", "s", "h", "hs")
.TYBITS <- c(1, 2, 4, 8, 16)
.TYLAST <- 32

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



# documentation of class sostpp -------------------------------------------

#'@name sostpp.object
#'@aliases sostpp.object sostpp
#'@title Class of Second-Order Stationarity-Typed Point Patterns
#'@description
#'  A class \code{"sostpp"} to represent a two-dimensional point
#'  pattern, with extra information about hidden second-order 
#'  stationarity. Built on top of the spatstat-class \code{"\link{ppp}"}.
#'@details
#'  This class represents a two-dimensional point pattern dataset with
#'  extra information relevant for analysis as a realization of a 
#'  second-order stationary point process.
#'  
#'  From its ancestor, spatstat-class \code{"\link{ppp}"}, an object of type 
#'  \code{sostpp} inherits the elements
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
#'    \cr\code{typemarks} \tab a data frame containing information relevant to 
#'    the type of stationarity,
#'    \cr\code{extra} \tab a container of additional information for internal use\cr
#'  }
#'  
#'  The type of second-order stationarity of an object \code{X} of class \code{sostpp}
#'  is returned as \code{character} by the function \code{\link{currenttype}}. Whether \code{X} has a given type of
#'  second-order stationarity, can be checked with \code{\link{has.type}}.
#'  
#'  Possible types of second-order stationarity (s.o.s.) are
#'  \tabular{ll}{
#'    \code{"w"}  \tab intensity reweighted s.o.s. 
#'                \cr\tab \code{marks} contains a row \code{lambda} with intensity evaluated
#'               in each data point
#'    \cr\code{"t"} \tab obtained by coordinate transformation
#'               \cr\tab \code{marks} contains rows \code{x0} and \code{y0}
#'                  of original (backtransformed) coordinates
#'                \cr\tab for each data point
#'    \cr\code{"s"} \tab locally rescaled s.o.s.
#'                  \cr\tab \code{marks} contains a row \code{invscale}, the inverse scale factor in each data point\cr
#'    \cr\code{"h"}  \tab homogeneous, to be evaluated with standard methods 
#'    \cr\code{"hs"} \tab homogeneous, to be evaluated with scale invariant statistics.  
#'    }
#'
#' 
#'@author Ute Hahn  \email{ute@@imf.au.dk}
#'
#'@keywords spatial attribute
#' 
NA 
