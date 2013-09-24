\name{sostpp.object}
\alias{sostpp.object}
\alias{sostpp}
\title{Class of Second-Order Stationarity-Typed Point Patterns}
\description{
  A class \code{"sostpp"} to represent a two-dimensional point
  pattern, with extra information about hidden second-order 
  stationarity. Built on top of the spatstat-class \code{"\link{ppp}"}.
}
\details{
  This class represents a two-dimensional point pattern dataset with
  extra information relevant for analysis as a realization of a 
  second-order stationary point process.
  
  From its ancestor, spatstat-class \code{"\link{ppp}"}, an object of type 
  \code{sostpp} inherits the elements
  \tabular{ll}{
    \code{x} \tab vector of \eqn{x} coordinates of data points \cr
    \code{y} \tab vector of \eqn{y} coordinates of data points \cr
    \code{n} \tab number of points \cr
    \code{window} \tab window of observation \cr
                  \tab (an object of class \code{\link{owin}}) \cr
    \code{marks} \tab vector or data frame of marks
  }
  Additionally, it contains the elements
  \tabular{ll}{
    \code{sostype} \tab integer, encrypts the type of stationarity that is assumed in analysis, see below \cr
    \code{typemarks} \tab a data frame containing information relevant to the type of stationarity,
    \cr
    \code{extra} \tab a container of additional information for internal use\cr
  }
  
  The type of second-order stationarity of an object \code{X} of class \code{sostpp}
  is returned as \code{character} by the function \code{\link{currenttype}}. Whether \code{X} has a given type of
  second-order stationarity, can be checked with \code{\link{has.type}}.
  
  Possible types of second-order stationarity (s.o.s.) are
  \tabular{ll}{
    \code{"w"}  \tab intensity reweighted s.o.s. \cr
               \tab \code{marks} contains a row \code{lambda} with intensity evaluated
               in each data point\cr
    \code{"t"} \tab obtained by coordinate transformation\cr
                  \tab \code{marks} contains rows \code{x0} and \code{y0}
                  of original (backtransformed) coordinates for each data point\cr
    \code{"s"} \tab locally rescaled s.o.s.\cr
                  \tab \code{marks} contains a row \code{invscale}, the inverse scale factor in each data point\cr
    \code{"h"}  \tab homogeneous, to be evaluated with standard methods \cr
    \code{"hs"} \tab homogeneous, to be evaluated with scale invariant statistics.  \cr
    }
}
 
\author{Ute Hahn  \email{ute@imf.au.dk}}

\keyword{spatial}
\keyword{attribute}
 
 
