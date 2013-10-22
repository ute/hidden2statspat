# plotting utilities and constants



# @param col color
# @param add coefficient on rgb, a single number
# @rdname sostatpp-internal
# @keywords internal
# @export
# @details works on all devices that support rgb

# .lightcoldefault <- .66

#' @title Plot of point pattern with test quadrats
#' @description Plot a point pattern with two sets of coloured test quadrats,
#' meant as a visualisation supplement for the function \code{\link{Kpermute.test}}.
#' @param x point pattern, object of class \code{\link{sostyppp}} or  class \code{spatstat::\link[spatstat]{ppp}}
#' @param quads1,quads2 test quadrats, a list of windows (\code{"\link{owin}"})
#' or a tessellation (\pkg{spatstat}-object \code{"\link[spatstat]{tess}"})
#' @param quads0 optional unused test quadrats
#' @param style1,style2,style0 plotting parameters, see details
#' @param bgcol1,bgcol2,bgcol0 colors, see details
#' @param light default lightness of the background, a number between 0 (dark) and 1 (white)
#' @param ... further arguments for plotting
#' @details One of the arguments \code{bgcol} or \code{style} has to be given.
#' If both are specified, \code{bgcol} overrides \code{style}. Lightness \code{light}
#' applies to the background colors, and is overridden by lightness specified in the \code{style}s.
# @param col color
# @param add coefficient on rgb, a single number
# @rdname sostatpp-internal
# @keywords internal
#' @export
#' @examples
#' # tile the beilschmiedia pattern evenly into many quadrats, and count points
#' data(bei)
#' allquads <- tiles(quadrats(bei, nx=8, ny=4))
#' npts <- sapply(allquads, function(w) npoints(bei[w]))
#'
#' # select quadrats with more than 30 points, and split into three sets
#' enoughpoints <- npts > 30
#' nmedian <- median(npts[enoughpoints])
#' few <- npts <= nmedian & enoughpoints
#' many <- npts >= nmedian
#'
#' quads1 <- allquads[many]
#' quads2 <- allquads[few]
#' quads0 <- allquads[!enoughpoints]
#'
#' redstyle <- list(col = "red", light = 0.5)
#' bluestyle <- list(col = "blue")
#' quadratsplot(bei, quads1, quads2, redstyle, bluestyle)
#'
#' # overriding some options
#' quadratsplot(bei, quads1, quads2, redstyle, bluestyle,
#'              light=.8, pch=16, cex=.4)
#'
#' # also unused quadrats
#' quadratsplot(bei, quads1, quads2, redstyle, bluestyle,
#'              light = .8, pch = 16, cex = .4,
#'              quads0 = quads0, bgcol0 = "yellow"
#'      )


oldquadratsplot <- function (x, quads1, quads2,
                          style1 = NULL, style2 = NULL,
                          bgcol1 = NULL, bgcol2 = NULL,
                          quads0 = NULL, style0 = list(col = "white"), bgcol0 = NULL,
                          light = 0.75,
                          ...)
{
  dotargs <- list(...)
  grapar <- par(no.readonly = TRUE)
  graphopt <- function (style, bgcol = NULL)
  {
    if (is.null(bgcol))
    {
      if (is.null(style)) stop(" at least one of bgcol and style has to be specified")
      sty <- uniquelist(c(dotargs, style))
      if (is.null(sty$col)) stop(" no color given ")
      bgcol <- if (is.null(sty$light)) lightcol(sty$col, light) else lightcol(sty$col, sty$light)
    }
    else bgcol <- lightcol(bgcol, light)
    grap <- updateoptions(updateoptions(grapar, style), dotargs)
    grap$col <- bgcol
    grap$add = TRUE
    return(grap)
  }

  grapho <- updateoptions(grapar, dotargs)
  # plot(x$window, main=grapho$main) this does not work
  if (is.null(dotargs$main)) plot(x$window, main = "") else plot(x$window, main = dotargs$main)
  if("tess" %in% class(quads1)) quads1 <- tiles(quads1)
  if("tess" %in% class(quads2)) quads2 <- tiles(quads2)
  if("tess" %in% class(quads0)) quads0 <- tiles(quads0)

  # don't plot dashed or dotted borders

  style1$lty <- "solid"
  style2$lty <- "solid"
  if (!is.null(style0)) style0$lty <- "solid"

  for (w in quads1) do.call(plot, c(list(w), graphopt(style1, bgcol1)))
  for (w in quads2) do.call(plot, c(list(w), graphopt(style2, bgcol2)))
  if (!is.null (quads0)) for (w in quads0) do.call(plot, c(list(w), graphopt(style0, bgcol0)))

  do.call(points, c(list(x), grapho))
}
