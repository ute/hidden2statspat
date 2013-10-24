#' Plot a second-order stationarity typed point pattern
#'
#' Extends the spatstat method \code{\link[spatstat]{plot.ppp}} by background color.
#'
#' @param x sos-typed point pattern which is plotted.
#' @param main title text
#' @param ... arguments passed to \code{\link{plot.ppp}}
#' @param col.win color of the plot window
#' @param alpha.win numeric between 0 and 1, alpha-value of window color, see
#'   \code{\link{alphacol}}
# @S3method plot sostyppp
#' @method plot sostyppp
#' @export plot.sostyppp
#' @seealso \code{\link{plot.ppp}} for the print method of class ancestor \code{ppp}.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

plot.sostyppp <- function(x, main, ..., col.win = NULL, alpha.win = 0.4) {
  dotargs <- style(...)
  unknownargs <- matching(dotargs, plot.owin, plot.ppp, plot.sostyppp, 
                          .graphparams, .plotparams, .notmatching = TRUE)
  if (length(unknownargs) > 0) 
    warning(paste("unused parameters: ", paste(names(unknownargs), collapse = ", ")))
  if (missing(main)) main <- deparse(substitute(x))
  allargs <- style(dotargs, main = main, col.win = col.win, alpha.win = alpha.win, 
                   NULL.rm = TRUE)
  if (!is.null(allargs$col.win)){
    # plot a coloured window first if plot is not added
      useargs <- matching(allargs, plot.owin, .graphparams, .plotparams) 
      useargs$col <- alphacol(allargs$col.win, allargs$alpha.win)
      do.call ("plot.owin", c(list(x$window), useargs))
      allargs$add <- TRUE
      allargs$main = NULL
  }
  # strip arguments that are only used for the window
  useargs <- matching(allargs, plot.ppp, .graphparams)
  useargs$col <- NULL
  do.call("plot.ppp", c(list(x), useargs))
}
