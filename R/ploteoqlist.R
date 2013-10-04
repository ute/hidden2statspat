# plot estonQuadrats

.plotparams <- list(type = "p",  xlim = NULL, ylim = NULL,
                    log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                    ann = par("ann"), axes = TRUE, frame.plot = TRUE,
                    panel.first = NULL, panel.last = NULL, asp = NA,
                    mgp = par("mgp"))

#' Plot a list of template K-functions estimated on quadrats
#'  
#' Plot a list of K-function estimates, as obtained by \code{\link{estOnQuadrats}}
#' @param x list of summary function estimates, object of class \code{eoqlist}
#' @param style optional list of plot options, as in \code{summaryplot.fdsample}
# @param col,colquad colors for plotting mean and single estimates
# @param lwd,lwdquad line widths
#' @param minn integer, minimal number of points on quadrat to allow for printing of
#'        the individual \eqn{K}-function, defaults to 10 
#' @param minm integer, minimal number of "valid" quadrats with at least \code{minn} points
# @param ylim y axis limits for plotting
#' @param ... further arguments for \code{\link{plot}},
#' @param add logical, whether to add to current plot. Defaults to FALSE.
# @param lwdtheo,ltytheo,coltheo plot parameters for K-function of Poisson point process (CSR).
#        If lwdtheo=0, the curve is not plotted. Defaults: lwdtheo = 1, ltythgeo="dotted", coltheo="black". 
#' @param theostyle plot options for the summary function of a Poisson point process (CSR), a list with
#' default elements \code{col}, \code{lwd} and \code{lty}, see \code{plot.fdsample}.
#' If \code{theostyle == NULL}, the CSR-curve is not plotted.
#' @param labline numerical, where to plot the axis labels
#' @details This method, as well as class \code{eoqlist}, is likely to be replaced in a future version
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @S3method plot eoqlist
#' @export

plot.eoqlist <- function(x,  
                         style = NULL,
                         corrections = "iso",       
                         minn = 10, minm = 1,
                         ..., add = FALSE, 
                         theostyle = list(lwd = 1, lty = "dotted", col = "black"),
                         labline = 2.4)
{
  OK <- x$npts >= minn
  m <- sum(OK)
  if (m < minm) stop(paste("not enough subpatterns with at least minn=", minn, "points"))
 # rr <- x$r    
  dotargs <- list(...)
  xopt <- getoptions(x$footheo, c(style, dotargs))
  allopt <- uniquelist(c(xopt, unusedoptions(xopt, dotargs)))
  plotargs <- updateoptions(.plotparams, allopt)
  
  if (is.null(corrections)) corrections <- names(x$fooli)
  ylim <- dotargs$ylim
  if (is.null(ylim))
  {
    yfo <- sapply(corrections, function(ckey) yrange(x$fooli[[ckey]]) )
    ylim <- range(c(yrange(x$footheo), range(yfo)))
  }
  textline <- par("mgp")
  textline[1] <- labline
  
  if (!add) {
    if (!is.null(theostyle)) { 
      plotargstheo <- uniquelist(c(theostyle, plotargs) )
      do.call(plot, c(list(x$footheo), list(plotargstheo), list(ylim = ylim, mgp = textline)));
      add <- TRUE }
  }
  
  
  # workaround to avoid lwd = 0
  if( is.null(dotargs$lwd) & is.null(style$lwd)) style$lwd <- par("lwd") 
  
  for (ckey in corrections)
    if (!is.null(x$fooli[[ckey]]))
    {
      summaryplot (x$fooli[[ckey]], ploptions = style, 
                   mgp = textline, ylim = ylim, ..., add = add)    
      add <- TRUE
    }  
}  

