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

# @param correction character
#' @rdname sosspp-internal
#' @keywords internal

correctionkey <- function (correction)
 pickoption("correction", correction,
                         c(none="none",
                           border="border",
                           "bord.modif"="bord.modif",
                           isotropic="iso",
                           Ripley="iso",
                           trans="trans",
                           translate="trans",
                           translation="trans",
                           best="iso"),
                         multi=TRUE) 


#' Template K-function, or Delta K_dir, estimated on quadrats
#'  
#' Estimate the template \eqn{K}-function, or the \eqn{Delta K}-function, 
#' on subsamples of a point pattern
#' 
#' @param X point pattern, object of  \bold{spatstat}-class \code{\link{ppp}}
#' @param type optional character, the type of second-order stationarity assumed: 
#'   \itemize{
#'        \item \code{"w"} reweighted
#'        \item \code{"t"} retransformed
#'        \item \code{"s"} locally rescaled
#'        \item \code{"h"} homogeneous, i.e. first order stationary
#'        \item \code{"hs"} homogeneous, but evaluated as scaled \eqn{K}-function
#'   }
#'   Only the first match is used. 
#' @param quads quadrats for subsampling \code{X}. A \code{list} of objects
#' of \bold{spatstat}-class \code{\link{owin}} or an object of \bold{spatstat}-class \code{\link{tess}}.
#' @param rmin optional, lower bound, defaults to 0,
#' @param rmax upper bound for the radius,
#' @param rlen optional, number of steps in argument vector, defaults to 256,
#' @param correction a character vector giving the edge correction type, may be
#'   one of \code{"border"},  \code{"isotropic"}, \code{"translate"}, \code{"none"}.
#' @param DeltaKdir logical, whether to return the Delta K_dir function instead 
#' of the template K-function  
#' @param ... further arguments for \code{\link{Khidden}} or \code{\link{DeltaKdir}}
#' @return An object of class \code{foolist}, which is a list with items
#' \itemize{
#'    \item \code{npts} number of points on the quadrats
#'    \item \code{r} arguments of the \eqn{K}-function
#'    \item \code{foomean} unweighted mean of the estimated \eqn{K}-functions
#'    \item \code{fooarray} array of estimates, of dimension \code{rlen} x number of quadrats
#'    \item \code{footheo} theoretical values under CSR (Poisson process)
#'    \item \code{type} the type character
#'    \item \code{correction} character, the correction used
#'    \item \code{xlab, ylab} labels for plotting
#' }
#' @seealso \code{\link{Khidden}}, \code{\link{DeltaKdir}}
#'    for estimation of the template \eqn{K}-function or of Delta K_dir.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @keywords robust
#' @keywords spatial
#' @keywords nonparametric
#' @keywords ts    

KestOnQuadrats <- function(X, type=NULL, quads,
                        rmin = 0, rmax = 1.25, rlen = 100,
                        correction = "isotropic",
                        DeltaKdir = FALSE,
                        ...)
{
  corr <- correction[1]
  ckey <- correctionkey(correction)
  if (is.null(type)) type <- getlasttype(X)
  if (length(type) == 0) type <- "w"
  tindex <-  which(.TYPES == type)
  
  pplist <- ppsplit(X, quads)
  ckey <- correctionkey(correction)
  rr <- seq(rmin, rmax, length.out=rlen)
  if (DeltaKdir){  
    Kar <- sapply(pplist, function(X) DeltaKdir(X, type = type, r = rr, correction = corr, ...)[[ckey]])
    ylab <- substitute(widehat(Delta*K)[dir]^(x)*(r), list(x=type))
    Ktheo <- rep(0, rlen)
  } 
  else {
    Kar <- sapply(pplist, function(X) Khidden(X, type = type, r = rr, correction = corr, ...)[[ckey]])
    ylab <- substitute(widehat(K)[0]^(x)*(r), list(x = type))
    Ktheo <- rr^2 * pi
  }  
  Kmean <- apply(Kar, 1, mean, na.rm=T)
  npts <- sapply(pplist, npoints)
  Klist <- list(npts = npts,
    r = rr,
    foomean = Kmean,
    fooarray = Kar,
    footheo = Ktheo,    
    type = type,
    correction = correction,
    ylab = ylab,
    xlab = "r"
    )
  attr(Klist, "class") <- "foolist"
  return(Klist)
}


# @param col color
# @param alpha mixing coefficient, alpha channel
#' @rdname sosspp-internal
#' @keywords internal

alphacol <- function(col, alpha = 0.4)  rgb(t(col2rgb(col)/256), alpha=alpha)

#' Plot a list of template K-functions estimated on quadrats
#'  
#' Plot a list of K-function estimates, as obtained by \code{\link{KhiddenOnQuadrats}}
#' @param klist list of \eqn{K}-function estimates
#' @param colmean,colquad colors for plotting mean and single estimates
#' @param lwdmean,lwdquad line widths
#' @param minn integer, minimal number of points on quadrat to allow for printing of
#'        the individual \eqn{K}-function, defaults to 10 
#' @param minm integer, minimal number of "valid" quadrats with at least \code{minn} points
#' @param ylim y axis limits for plotting
#' @param ... further arguments for \code{\link{plot}},
#' @param add logical, whether to add to current plot. Defaults to FALSE.
#' @param lwdtheo,ltytheo plot parameters for K-function of Poisson point process (CSR).
#'        If lwdtheo=0, the curve is not plotted. Defaults: lwdtheo = 1, ltythgeo="dashed". 
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @S3method plot foolist
#' @export

plot.foolist <- function(flist,  
                  col = "red", colquad = alphacol(col, 0.4),
                  lwd = 2, lwdquad = 0.7, 
                  minn = 10, minm = 1,
                  ylim = NULL,
                  ..., add = FALSE, 
                  lwdtheo = 1, ltytheo = "dashed")
{
  OK <- flist$npts >= minn
  m <- sum(OK)
  if (m < minm) stop(paste("not enough subpatterns with at least minn=", minn, "points"))
  foomean <- apply(flist$fooarray[, OK], 1, mean, na.rm=T)  
  rr <- flist$r    
  if (is.null(ylim)) ylim <- c(min(flist$footheo, flist$fpparray), max(flist$fooarray, flist$footheo))

  if (!add) plot(rr, foomean, type = "l", ylim = ylim,
                xlab = flist$xlab, ylab = flist$ylab,
                col = col, lwd = lwd, ...)
  else lines(rr, foomean, col = col, lwd = lwd, ...)

  # plot Poisson reference curve
  if (!add & lwdtheo>0) lines(rr, flist$footheo, lwd = lwdtheo, lty = ltytheo)
  
  if(lwdquad > 0) apply(flist$fooarray, 2, function(y) lines(rr, y, col=colquad, lwd=lwdquad, ...))
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
#' @keywords spatial
#' @keywords nonparametric
#' @keywords ts    

Kpermtest <- function(X, Y = NULL,
                      type = NULL,
                      quad1, quad2,
                      rmin = 0,
                      rmax,
                      rlen = 100,
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
  ckey <- correctionkey(correction)
  rr <- seq(rmin, rmax, length.out=rlen)
  Kar1 <- sapply(pplist1, function(X) Khidden(X, type = type, r = rr, correction = corr, ...)[[ckey]])
  Kar2 <- sapply(pplist2, function(X) Khidden(X, type = type, r = rr, correction = corr, ...)[[ckey]])
  if (use.tbar) testerg <- studpermut.test((Kar1 / rr)[rr>0, ], (Kar2 / rr)[rr>0, ], use.tbar = TRUE, nperm = nperm)
  else  testerg <- studpermut.test(Kar1, Kar2, use.tbar = FALSE, nperm = nperm)
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