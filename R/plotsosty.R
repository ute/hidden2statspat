#' Plot a second-order stationarity typed point pattern
#'
#' Extends spatstat method \code{\link{plot.ppp}} by background color.
#'
#' @param x sos-typed point pattern which is plotted.
#' @param main title text 
#' @param ... arguments passed to \code{\link{plot.ppp}}
#' @param col.win color of the plot window
#' @param alpha.win numeric between 0 and 1, alpha-value of window color, see
#'   \code{\link{alphacol}}
#' @S3method plot sostyppp
#' @method plot sostyppp
# @export
#' @seealso \code{\link{plot.ppp}} for the print method of class ancestor \code{ppp}.
#' @author Ute Hahn,  \email{ute@@imf.au.dk}

plot.sostyppp <- function(x, main, ..., col.win = NULL, alpha.win = 0.4) {
  dotargs <- style(...)
  if (missing(main)) main <- deparse(substitute(x))
  if (!is.null(col.win)){
    # plot a coloured window first if plot is not added
      splot(x$window, main, col = alphacol(col.win, alpha.win), dotargs)
      # remove plot.owin specific parameter
      dotargs <- matching(dotargs, plot.owin, .notmatching = TRUE)
      class(dotargs) <- c("style", "list")
      dotargs$add <- TRUE
      main = NULL
  }
  splot(as.ppp(x), main, dotargs)
}