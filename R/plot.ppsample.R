#' Plot a sample of point patterns
#'
#' Plots all point patterns contained in an \code{ppsample} object.
#'
#' @param x the \code{ppsample} object that is to be plotted.
#' @param ... arguments passed to the plot methods of the point patterns in \code{x}.
#' @param add logical, if \code{TRUE} add to the current plot (default: \code{FALSE})
#' @param allinone logical, if \code{FALSE} start a new plot for every pattern (default: \code{TRUE}).
# @param allwindows logical, if \code{TRUE} plot the window for every pattern (default: \code{TRUE}).
# so ein quatsch, das ist dochn sample, mann! Also immer mit Rand
#' @S3method plot ppsample
#' @method plot ppsample
#' @export
#' @seealso \code{\link{plot.ppp}} and \code{\link{plot.sostyppp}}
#' for the plot methods of elements in the sample.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

plot.ppsample <- function(x, ..., add = FALSE, allinone = TRUE) {
  dotargs <- simplist(...)
  if (!add){
    # spatstatpatch
    # don't let plot.ppp or plot.owin print a rubbish name as title because of do.call
    if (is.null(dotargs$main)) {
      dotargs$main = ""
 #     dotargs$main <- deparse(substitute(x))[1] also gives rubbish,
  #    if plot.ppsample was invoked by do call
  #
 #### TODO: simplist$ wenn mal ganz viel Zeit ist... #################
      # class(dotargs) <- c("class", "list")
    }
    if (allinone) {
      w <- attr(x, "parentwindow")
      splot(w, dotargs)
    }
  }
  add <- allinone
  # force plotting the windows, even if no background colour given
  if (is.null(dotargs$col.win))
    dotargs$col.win <- par("bg")
  # remove eventual pure col arguments, since plot.ppp plots points in same color as window
  dotargs$col <- NULL
  do.call(lplot, c(list(x), dotargs, add=add, allinone = allinone))
}


#'splot method for class ppsample
#'
#'Plots an object of class \code{ppsample}. The function \code{splot.ppsample}
#'is an alias for \code{\link{plot.ppsample}}
#' @param x the \code{ppsample} object that is to be plotted.
#' @param ... arguments passed to the plot methods of the point patterns in \code{x}.
#' @param add logical, if \code{TRUE} add to the current plot (default: \code{FALSE})
#' @param allinone logical, if \code{FALSE} start a new plot for every pattern (default: \code{TRUE}).
#' @S3method splot ppsample
#' @method splot ppsample
#' @seealso \code{\link{plot.ppsample}}.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

splot.ppsample <- function(x, ..., add = FALSE, allinone = TRUE) {
  plot.ppsample(x, ..., add = add, allinone = allinone)
}