# Testing based on K-function on quadrats

# @param X point pattern
# @param quadrats a list of owins, or an object of class \code{"\link{tess}"}
#' @rdname sostatpp-internal
#' @keywords internal

ppsplit <- function (X, quadrats)
{
  # names throw warnings in plyr. Remove all names that are not list names from quads  
  unnamelist <- function(l)
{
  if (!is.list(l)) names(l) <- NULL
  if (is.list(l)) for (i in 1:length(l)) l[[i]] <- unnamelist(l[[i]])
  return(l)
}

  if("tess" %in% class(quadrats)) quadrats <- tiles(quadrats)
  nonamequads <- unnamelist(quadrats)
  return(lapply(nonamequads, function(w) X[w]))
}

# @param correction character
#' @rdname sostatpp-internal
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


#' Estimate summary function on quadrats
#'  
#' Estimate  a  summary function on subsamples of a point pattern.
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
# In the current version not used for retransformed used if \code{type != "t"}.
# @param tquads used instead of \code{quads} for the backtransformed pattern if \code{type != "t"}.
#' @param fun the summary function to be applied 
#' @param rmin optional, lower bound, defaults to 0,
#' @param rmax upper bound for the radius,
#' @param rlen optional, number of steps in argument vector, defaults to 256,
#' @param ... further arguments passed to \code{"fun"} or for type-tagging of the point pattern \code{X}
#' @return An object of class \code{foolist}, which is a list with items
#' \itemize{
#'    \item \code{npts} number of points on the quadrats
#'    \item \code{r} arguments of the \eqn{K}-function
#    \item \code{foomean} unweighted mean of the estimated \eqn{K}-functions
#'    \item \code{fooarray} array of estimates, of dimension \code{rlen} x number of quadrats
#'    \item \code{footheo} theoretical values under CSR (Poisson process)
#'    \item \code{type} the type character
#'    \item \code{correction} character, the correction used
#'    \item \code{xlab, ylab} labels for plotting
#' }
#' @seealso \code{\link{estK}}, \code{\link{estDeltaKdir}}
#'    for estimation of the template \eqn{K}-function or of Delta K_dir.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @keywords robust
#' @keywords spatial
#' @keywords nonparametric
#' @keywords ts    

estOnQuadrats <- function(X, type = NULL, quads, 
                          fun = estK,
                          rmin = 0, rmax = 1.25, rlen = 100,
                        ...)
{
  # corr <- correction[1]
  # ckey <- correctionkey(correction)
  
  if (!is.null(type)) # type request
  {
    type <- type[1]
    if (!has.type(X, type)) X <- as.sostyppp(X, type, ...)
  }
  else type <- currenttype(X)
  
  # now we are sure that X has the right type :-)
  
  if (length(type) == 0) stop ("no type of hidden 2nd order stationarity given")
  
  tindex <-  which(.TYPES == type)
  
  # names throw warnings in plyr. Remove all names from quads
  pplist <- ppsplit(X, quads)
  npp <- length(pplist)
  rr <- seq(rmin, rmax, length.out=rlen)

  # KK serves to get relevant information from fv, a tiny waste of time...
  KK <- fun(pplist[[1]], r = rr,  ...)
  
  xlab <- attr(KK, "xlab")
  ylab <- attr(KK, "ylab")
  if (is.call(ylab)) ylab <- do.call(expression, list(ylab))
  KKu <- attr(KK, "units")$singular
  xlab <- paste(attr(KK, "argu"), ifelse(KKu != "unit", paste("(", KKu, ")", sep=""), ""), sep=" ")
  dotnames <- attr(KK, "dotnames")
  estnames <- dotnames[! (dotnames %in% c(".id", "r", "theo"))]
  
  #Kar <- sapply(pplist, function(X) Kfun(X, r = rr, correction = corr, ...)[[ckey]])
  # have read that sapply is inefficient (?)
  
  Kar <- ldply(pplist, function(X) fun(X, r = rr, ...))
  Kars <- llply(estnames, function(x) 
    fdsample(rr, array(Kar[[x]], c(rlen, npp)), xlab = xlab, ylab=ylab))             
  names(Kars) <- estnames

  # CSR-function and plot axis labels from info in KK
  Ktheo <- fdsample(rr, KK$theo, xlab = xlab, ylab=ylab)
  
  npts <- sapply(pplist, npoints)
               
# <<< change  the structure later !!! >>> 
  Klist <- list(npts = npts,
   # r = rr,
  #  foomean = Kmean,
    # fooarray = Kar,
    fooli = Kars,            
    footheo = Ktheo,    
    type = type
 #   ylab = ylab,
  #  xlab = xlab
    )
  attr(Klist, "class") <- "eoqlist"
  return(Klist)
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
#' @param quads1,quads2 quadrats for subsampling \code{X} (and \code{Y}). A \code{list} of objects
#' of \bold{spatstat}-class \code{\link{owin}} or an object of \bold{spatstat}-class \code{\link{tess}}.
# In the current version not used for retransformed used if \code{type != "t"}.
# @param tquads1,tquads2 used instead of \code{quads1} and \code{quads2} for 
#      the backtransformed pattern if \code{type != "t"}.
#' @param rmin optional, lower integration bound, defaults to 0; see `Details',
#' @param rmax upper integration bound, see `Details',
#' @param rlen optional, number of steps for numerical integration, defaults to 256; see `Details',
#' @param Kfun optional \code{function}, the \eqn{K}-function to be used, 
#'  either \code{\link{estK}} (default) or \code{\link{estDeltaKdir}}
#' @param correction a character vector giving the edge correction type, may be
#'   any subset of \code{"border"},  \code{"isotropic"}, \code{"translate"}, \code{"none"}.
#' @param ... further arguments for \code{Kfun}
#' @param use.tbar logical, defaults to \code{FALSE}. Whether to apply 
#'   studentization after integration, see `Details'.
#' @param nperm \code{NULL} or an integer giving the number of random 
#'   permutations. If \code{NULL}, all permutations are used. Only feasible if 
#'   group sizes \eqn{m_1}, \eqn{m_2} do not exceed 10.
#' @details 
#'   if \code{Y} is given, the \eqn{K}-functions of the two point patterns \code{X} and {Y} 
#'   are compared, otherwise a test of hidden second-order stationarity is carried out on \code{X}.
#'   
#'   The test is a permutation test, as proposed in Hahn(2012) and Hahn & Jensen (2013). 
#'   Which test statistic is used, is controlled by 
#'   \code{use.tbar}:
#'     \itemize{
#'   \item \eqn{T = mean [ (K_1(r)-K_2(r))^2 / (s_1^2(r)/m_1 + s_2^2(r)/m_2), ]} if \code{use.tbar==FALSE}    or
#'   \item \eqn{Tbar = mean [ (K_1(r)-K_2(r))^2 ] / mean[ s_1^2(r)/m_1 + s_2^2(r)/m_2, ]} if\code{use.tbar==TRUE},
#'   }
#'   where \eqn{K_1(r)} and \eqn{K_2(r)} are the group means, and 
#'    \eqn{s_1^2(r), s_2^2(r)} are within group 
#'   variances at a fixed argument \eqn{r}. 
#'   The mean is taken over the interval \code{[rmin, rmax]}.
#'   Note that test statistic \eqn{T} differs slightly from the statistics
#'   given in the above papers, in that the integral is replaced by the mean here.
#'   
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
# @source Hahn(2012), with slight modification (using the mean instead of the 
#    integral, in order to avoid having to pass the arguments of the functions)
#' @references Hahn, U. (2012) A Studentized Permutation Test for the Comparison 
#' of Spatial Point Patterns. \emph{Journal of the American Statistical 
#' Association},  \bold{107} (498), 754--764.
#' @references Hahn, U. and Jensen, E. B. V. (2013)
#' Inhomogeneous spatial point processes with hidden second-order stationarity.
#' \emph{CSGB preprint} 2013-07, \url{http://data.imf.au.dk/publications/csgb/2013/math-csgb-2013-07.pdf}
#' @keywords htest
#' @keywords robust
#' @keywords spatial
#' @keywords nonparametric
#' @keywords ts    

Kpermute.test <- function(X, Y = NULL,
                      type = NULL,
                      quads1, quads2,
                    #  tquads1, tquads2,
                      rmin = 0,
                      rmax,
                      rlen = 100,
                      Kfun = estK, 
                      correction = "iso",
                      ...,
                      use.tbar = FALSE,
                      nperm = 25000)
{
  AnisTest <- identical(Kfun, estDeltaKdir)
  
  dataname <- if(is.null(Y)) paste( "point pattern",deparse(substitute(X)))
              else paste( "point patterns",deparse(substitute(X)), "and" ,deparse(substitute(Y)))
    
  if (!is.null(type)) # type request
  {
    type <- type[1]
    if (!has.type(X, type)) X <- as.sostyppp(X, type, ...)
  }
  else type <- currenttype(X)
  
  # now we are sure that X has the right type :-)
  
  if (length(type) == 0) stop ("no type of hidden 2nd order stationarity given")
  
  tindex <-  which(.TYPES == type)
  
  pplist1 <- ppsplit(X, quads1)
  pplist2 <-  if(is.null(Y)) ppsplit(X, quads2) else ppsplit(Y, quads2)
  corr <- correction[1]
  ckey <- correctionkey(correction)
  rr <- seq(rmin, rmax, length.out=rlen)
  Kar1 <- sapply(pplist1, function(X) Kfun(X, r = rr, correction = corr, ...)[[ckey]])
  Kar2 <- sapply(pplist2, function(X) Kfun(X, r = rr, correction = corr, ...)[[ckey]])
  if (use.tbar)
  {  
    Ksample1 <- fdsample(rr[rr>0], (Kar1 / rr)[rr>0, ])
    Ksample2 <- fdsample(rr[rr>0], (Kar2 / rr)[rr>0, ])
  } 
  else
  {  
    Ksample1 <- fdsample(rr, Kar1)
    Ksample2 <- fdsample(rr, Kar2)
  } 
  testerg <-tL2.permtest(Ksample1, Ksample2, nperm = nperm, use.tbar = use.tbar) 
  if (is.null(Y))  {
    method <- c(paste("Studentized permutation test of",
                              .TYPENAMESX[tindex]," hidden second-order stationarity,"),
                        ifelse(AnisTest, "directional version, using Delta K_dir",
                               "isotropic version, using K_0"),
                        paste("test statistic: ", if(use.tbar) "Tbar," 
                              else "T,", "integration limits",rmin,"to",rmax),
                        testerg$method[2] )
    alternative <- c(paste("not the same",.TYPENAMESX[tindex],"K-function"),
                                 if(AnisTest) ",\nbut different kinds of anisotropy")
  } else  {
    method <- c(paste("Studentized permutation test of identical",
                              .TYPENAMESX[tindex]," K-functions,"),
                        ifelse(AnisTest, "directional version, using Delta K_dir",
                               "isotropic version, using K_0"),
                        paste("test statistic: ", if(use.tbar) "Tbar," 
                              else "T,", "integration limits",rmin,"to",rmax),
                        #if(is.null(nperm)) "exact test, using all permutations" 
                        #else paste("using",nperm,"randomly selected permuations")
                        testerg$method[2])
    alternative <- c(paste("not the same",.TYPENAMESX[tindex],"K-function"),
                                 if(AnisTest) ",\nbut different kinds of anisotropy")
  }
  testerg$method <- method
  testerg$alternative <- alternative
  testerg$data.name <- dataname
  return(testerg)
}

