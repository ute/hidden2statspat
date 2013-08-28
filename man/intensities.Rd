\name{intensities}
\alias{intensities}
\alias{bronze.intens}
\alias{bei.intens.bayes}
\alias{bei.intens.maxlik}
\alias{bei.intens.nonpar}
\alias{scholtzia.intens}
\docType{data}
\title{Fitted intensities}
\description{
  Intensities of the inhomogeneous point patterns \code{\link{bei}}%, 
  %\code{\link{bronzefilter}} 
  and \code{\link{scholtzia}}, as used in the
  paper Hahn and Jensen (2013). These intensities are included in the package 
  for reason of reference.
} 
\format{
    
  \code{bei.intens.bayes}, \code{bei.intens.maxlik} and \code{bei.intens.nonpar}
  are 101x201 pixel images, representing Bayesian, Maximum Likelihood and a nonparametric
  estimate of the intensity of the point pattern \code{\link{bei}}.
  See \code{\link{im.object}} for details of the format.
  
  \code{bronze.intens} is a  128x128 pixel image representing the estimated 
  intensity of the point pattern \code{\link{bronzefilter}}. 
  
  
  
  \code{scholtzia.intens} is a 128x128 pixel image representing the estimated 
  intensity of the point pattern \code{\link{scholtzia}}.
  
}
\usage{data(intensities)}
\source{
  Hahn and Jensen (2013).
}
\references{
  Hahn, U. and Jensen, E.B.V. (2013)
  Inhomogeneous spatial point processes with hidden second-order stationarity.
} 
\keyword{datasets}
\keyword{spatial}

