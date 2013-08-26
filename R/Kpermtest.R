# Testing based on K-function on quadrats

# @param X point pattern
# @param quadrats a list of owins, or an object of class \code{"\link{tess}"}
#' @rdname sosspp-internal
#' @keywords internal

ppsplit <- function (X, quadrats)
{
  if("tess" %in% class(quadrats)) quadrats <- tiles(quadrats)
  return(lapply(quadrats, function(w) X[w]))
}


# from MakeHidden:
.TYPENAMES  <- c("reweighted", "retransformed", "rescaled", "homogeneous", "homogeneous, scaled")
.TYPENAMESX  <- c("reWeighted", "reTransformed", "reScaled", "Homogeneous", "Homogeneous, Scaled")
.TYPES  <- c("w", "t", "s", "h", "hs")
.TYBITS <- c(1, 2, 4, 8, 16)
.TYLAST <- 32

#' Studentized Permutation Test for K-functions of Subsampled Point Patterns
#' 
#' Perform a studentized permutation test of equal \eqn{K}-function for one or two 
#' subsampled point patterns.
#' @param X point pattern, object of  \bold{spatstat}-class \code{\link{ppp}}
#' @param Y optional; point pattern, object of  \bold{spatstat}-class \code{\link{ppp}}
#' @param type optional character, the type of second-order stationarity assumed: 
#'   \itemize{
#'        \item \code{"w"} reweighted
#'        \item \code{"t"} retransformed
#'        \item \code{"s"} locally rescaled
#'        \item \code{"h"} homogeneous, i.e. first order stationary
#'        \item \code{"hs"} homogeneous, but evaluated as scaled \eqn{K}-function
#'   }
#'   Only the first match is used. 
#' @param quad1,quad2 quadrats for subsampling \code{X} (and \code{Y}). A \code{list} of objects
#' of \bold{spatstat}-class \code{\link{owin}} or an object of \bold{spatstat}-class \code{\link{tess}}.
#' @param rmin optional, lower integration bound, defaults to 0; see `Details',
#' @param rmax upper integration bound, see `Details',
#' @param rlen optional, number of steps for numerical integration, defaults to 256; see `Details',
#' @param correction a character vector giving the edge correction type, may be
#'   any subset of \code{"border"},  \code{"isotropic"}, \code{"translate"}, \code{"none"}.
#' @param ... further arguments for \code{\link{Khidden}}
#' @param use.tbar logical, defaults to \code{FALSE}. Whether to apply 
#'   studentization after integration, see `Details'.
#' @param nperm \code{NULL} or an integer giving the number of random 
#'   permutations. If \code{NULL}, all permutations are used. Only feasible if 
#'   group sizes \eqn{m_1}, \eqn{m_2} do not exceed 10.
#' @details 
#'   if \code{Y} is given, the \eqn{K}-functions of the two point patterns \code{X} and {Y} 
#'   are compared, otherwise a test of hidden second-order stationarity is carried out on \code{X}.
#'   
#'   !!!!!!!!!!!!!!!!!!!!! change this !!!!!!!!
#'   Test statistics are integral means of studentized square distances 
#'   between the group means. The test statistics are closely related, but not 
#'   equal to Hotelling's two-sample T-squared statistic. It is assumed that 
#'   the functions all have the same, equidistant arguments. Depending on the 
#'   value of \code{use.tbar}, the test statistic is either
#'   \itemize{
#'   \item \eqn{T = mean [ (mu_1(x)-mu_2(x))^2 / (s_1^2(x)/m_1 + s_2^2(x)/m_2) ]}    or
#'   \item \eqn{Tbar = mean [ (mu_1(x)-mu_2(x))^2 ] / mean[ s_1^2(x)/m_1 + s_2^2(x)/m_2 ]}
#'   }
#'   where \eqn{mu_1(x), mu_2(x)} and \eqn{s_1^2(x), s_2^2(x)} are within group 
#'   means and variances at a fixed argument \eqn{x}, and the overall
#'   the mean is taken over all \eqn{n} arguments \eqn{x}.
#'   
#'   If \code{nperm == NULL}, the exact test with all permutations is used
#'   (combinations, for symmetry reasons). This may cause memory or computing 
#'   time issues.
#'   If \code{nperm} is given as an integer, the permutations are \code{sample}d randomly, 
#'   unless \code{nperm} is larger than the number of disjoint combinations. In that
#'   case, the exact version is used. 
#' @return 
#' A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the test statistic,}
#' \item{p.value}{the p-value of the test,}
#' \item{alternative}{a character string describing the alternative hypothesis,}
#' \item{method}{a character string indicating what type of test was performed,}
#' \item{data.name}{a character string giving the name(s) of the data.}
#'    
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @source Hahn(2012), with slight modification (using the mean instead of the 
#'    integral, in order to avoid having to pass the arguments of the functions)
#' @references Hahn, U. (2012) A Studentized Permutation Test for the Comparison 
#' of Spatial Point Patterns. \emph{Journal of the American Statistical 
#' Association},  \bold{107} (498), 754--764.
#' @references Hahn, U. and Jensen, E. B. V. (2013)
#' Inhomogeneous spatial point processes with hidden second-order stationarity.
#' @keywords htest
#' @keywords robust
#' @keywords nonparametric
#' @keywords ts    

Kpermtest <- function(X, Y = NULL,
                      type = NULL,
                      quad1, quad2,
                      rmin = 0,
                      rmax,
                      rlen = 256,
                      correction = "iso",
                      ...,
                      use.tbar = FALSE,
                      nperm = NULL)
{
  if (is.null(type)) type <- getlasttype(X)
  if (length(type) == 0) type <- "w"
  tindex <-  which(.TYPES == type)
  
  pplist1 <- ppsplit(X, quads1)
  pplist2 <-  if(is.null(Y)) ppsplit(X, quads2) else ppsplit(Y, quads2)
  corr <- correction[1]
  rr <- seq(rmin, rmax, length.out=rlen)
  Kar1 <- sapply(pplist1, function(X) Khidden(X, type = type, r = rr, correction = corr, ...)[[corr]])
  Kar2 <- sapply(pplist2, function(X) Khidden(X, type = type, r = rr, correction = corr, ...)[[corr]])
  testerg <- studpermut.test(Kar1, Kar2, use.tbar = use.tbar, nperm = nperm)
  if (is.null(Y))  {
    testerg$method <- c(paste("Studentized permutation test of",
                              .TYPENAMESX[tindex]," hidden second-order stationarity,"),
                        paste("test statistic: ", if(use.tbar) "Tbar," 
                              else "T,", "integration limits",rmin,"to",rmax),
                        if(is.null(nperm))"exact test, using all permutations" 
                        else paste("using",nperm,"randomly selected permuations") )
    testerg$alternative <- paste("not",.TYPENAMESX[tindex],"second-order stationary")
    testerg$data.name <-  paste( "point pattern",deparse(substitute(X)))
  } else  {
    testerg$method <- c(paste("Studentized permutation test of identical",
                              .TYPENAMESX[tindex]," K-functions,"),
                        paste("test statistic: ", if(use.tbar) "Tbar," 
                              else "T,", "integration limits",rmin,"to",rmax),
                        if(is.null(nperm)) "exact test, using all permutations" 
                        else paste("using",nperm,"randomly selected permuations"))
    testerg$alternative <- paste("not the same",.TYPENAMESX[tindex],"K-function")
    testerg$data.name <-  paste( "point patterns",deparse(substitute(X)), "and" ,deparse(substitute(Y)))
  }
  return(testerg)
}



#' Studentized Permutation Test of Stationarity in Subsampled Point Patterns
#' 
#' Perform a studentized permutation test of equal \eqn{K}-function for one or two 
#' subsampled point patterns.
#' @param X point pattern, object of  \bold{spatstat}-class \code{\link{ppp}}
# @param Y optional; point pattern, object of  \bold{spatstat}-class \code{\link{ppp}}
# @param type optional character, the type of second-order stationarity assumed: 
#   \itemize{
#        \item \code{"w"} reweighted
#        \item \code{"t"} retransformed
#        \item \code{"s"} locally rescaled
#        \item \code{"h"} homogeneous, i.e. first order stationary
#        \item \code{"hs"} homogeneous, but evaluated as scaled \eqn{K}-function
#   }
#   Only the first match is used. 
#' @param quad1,quad2 quadrats for subsampling \code{X}. A \code{list} of objects
#' of \bold{spatstat}-class \code{\link{owin}} or an object of \bold{spatstat}-class \code{\link{tess}}.
#' @param rmin optional, lower integration bound, defaults to 0; see `Details',
#' @param rmax upper integration bound, see `Details',
#' @param rlen optional, number of steps for numerical integration, defaults to 256; see `Details',
# @param correction a character vector giving the edge correction type, may be
#   any subset of \code{"border"},  \code{"isotropic"}, \code{"translate"}, \code{"none"}.
#' @param ... further arguments for \code{\link{DeltaKdir}}
#' @param use.tbar logical, defaults to \code{FALSE}. Whether to apply 
#'   studentization after integration, see `Details'.
#' @param nperm \code{NULL} or an integer giving the number of random 
#'   permutations. If \code{NULL}, all permutations are used. Only feasible if 
#'   group sizes \eqn{m_1}, \eqn{m_2} do not exceed 10.
#' @details 
#   if \code{Y} is given, the \eqn{K}-functions of the two point patterns \code{X} and {Y} 
#   are compared, otherwise a test of hidden second-order stationarity is carried out on \code{X}.
#   
#'   !!!!!!!!!!!!!!!!!!!!! change this !!!!!!!!
#'   Test statistics are integral means of studentized square distances 
#'   between the group means. The test statistics are closely related, but not 
#'   equal to Hotelling's two-sample T-squared statistic. It is assumed that 
#'   the functions all have the same, equidistant arguments. Depending on the 
#'   value of \code{use.tbar}, the test statistic is either
#'   \itemize{
#'   \item \eqn{T = mean [ (mu_1(x)-mu_2(x))^2 / (s_1^2(x)/m_1 + s_2^2(x)/m_2) ]}    or
#'   \item \eqn{Tbar = mean [ (mu_1(x)-mu_2(x))^2 ] / mean[ s_1^2(x)/m_1 + s_2^2(x)/m_2 ]}
#'   }
#'   where \eqn{mu_1(x), mu_2(x)} and \eqn{s_1^2(x), s_2^2(x)} are within group 
#'   means and variances at a fixed argument \eqn{x}, and the overall
#'   the mean is taken over all \eqn{n} arguments \eqn{x}.
#'   
#'   If \code{nperm == NULL}, the exact test with all permutations is used
#'   (combinations, for symmetry reasons). This may cause memory or computing 
#'   time issues.
#'   If \code{nperm} is given as an integer, the permutations are \code{sample}d randomly, 
#'   unless \code{nperm} is larger than the number of disjoint combinations. In that
#'   case, the exact version is used. 
#' @return 
#' A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the test statistic,}
#' \item{p.value}{the p-value of the test,}
#' \item{alternative}{a character string describing the alternative hypothesis,}
#' \item{method}{a character string indicating what type of test was performed,}
#' \item{data.name}{a character string giving the name(s) of the data.}
#'    
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @source Hahn and Jensen (2013)
#' @references Hahn, U. and Jensen, E. B. V. (2013)
#' Inhomogeneous spatial point processes with hidden second-order stationarity.
#' @keywords htest
#' @keywords robust
#' @keywords nonparametric
#' @keywords ts    

Kanistest <- function(X, #Y = NULL,
                      # type = NULL,
                      quad1, quad2,
                      rmin = 0,
                      rmax,
                      rlen = 256,
                      # correction = "iso",
                      ...,
                      use.tbar = FALSE,
                      nperm = NULL)
{
  # if (is.null(type)) type <- getlasttype(X)
  type <- "s"
  tindex <-  which(.TYPES == type)
  
  pplist1 <- ppsplit(X, quads1)
  pplist2 <-  if(is.null(Y)) ppsplit(X, quads2) else ppsplit(Y, quads2)
  corr <- "iso" #correction[1]
  rr <- seq(rmin, rmax, length.out=rlen)
  Kar1 <- sapply(pplist1, function(X) DeltaKdir(X, type = type, r = rr, correction = corr, ...)[[corr]])
  Kar2 <- sapply(pplist2, function(X) DeltaKdir(X, type = type, r = rr, correction = corr, ...)[[corr]])
  testerg <- studpermut.test(Kar1, Kar2, use.tbar = use.tbar, nperm = nperm)
  #if (is.null(Y))
    {
    testerg$method <- c(paste("Studentized permutation test of stationarity, based on scaled K-functions"),
                        paste("test statistic: ", if(use.tbar) "Tbar,"  else "T,", 
                              "integration limits",rmin,"to",rmax),
                        if(is.null(nperm))"exact test, using all permutations" 
                        else paste("using",nperm,"randomly selected permuations") )
    testerg$alternative <- paste("not stationary, but different kinds of anisotropy")
    testerg$data.name <-  paste( "point pattern",deparse(substitute(X)))
  }# else  {
  #  testerg$method <- c(paste("Studentized permutation test of identical",
  #                            .TYPENAMESX[tindex]," K-functions,"),
  #                      paste("test statistic: ", if(use.tbar) "Tbar," 
  #                            else "T,", "integration limits",rmin,"to",rmax),
  #                      if(is.null(nperm)) "exact test, using all permutations" 
  #                      else paste("using",nperm,"randomly selected permuations"))
  #  testerg$alternative <- paste("not the same",.TYPENAMESX[tindex],"K-function")
  #  testerg$data.name <-  paste( "point patterns",deparse(substitute(X)), "and" ,deparse(substitute(Y)))
  #}
  return(testerg)
}