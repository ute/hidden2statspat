# plotting utilities and constants



# @param col color
# @param add coefficient on rgb, a single number
# @rdname sostatpp-internal
# @keywords internal
# @export
# @details works on all devices that support rgb

# .lightcoldefault <- .66

# @param col color
# @param add coefficient on rgb, a single number
# @rdname sostatpp-internal
# @keywords internal
# @export
# @details works on all devices that support rgb

# lightcol <- function(col, light = .lightcoldefault){
#  RGB <- pmin(light + (1 - light) * col2rgb(col) / 255, 1) # pmin just to be on the safe side
#  rgb(RGB["red", ], RGB["green", ], RGB["blue", ])
#}  

# @param col color
# @param add coefficient on rgb, a single number
#' @rdname sostatpp-internal
#' @keywords internal
#' @export
# @details works on all devices that support rgb

quadratsplot <- function (X, quads1, quads2, col1, col2, main="", ...)
{
  plot(X$window, main = main, ...)
  if("tess" %in% class(quads1)) quads1 <- tiles(quads1)
  if("tess" %in% class(quads2)) quads2 <- tiles(quads2)
  for (w in quads1) plot(w, col=col1, add=T, ...)
  for (w in quads2) plot(w, col=col2, add=T, ...)
  points(X, ...)
}
