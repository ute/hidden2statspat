#' @useDynLib sosspp
# @param x,y vectors of x and y coordinates
# @param lambda intensity vector
# @param invscale inverse scale factor
#' @rdname sosspp-internal
#' @keywords internal
#' @export
#'


lspairdist <- function(x, y, invscale, lambda)
{
  if (missing(invscale)) .Call("lspairdistla", x, y, lambda)
  else .Call("lspairdist", x, y, invscale)
}
  

# @param X point pattern
# @param ordered logical, is ignored
# @param invscale a vector of inverse scale factors in the points of X
# @param lambda a vector of intensities in the points of X
# @param rmax maximum locally scaled distance 
# @param ordered logical, is ignored
# @returns the locally scaled close pairs, where d= Euclidean distance, 
#   lsd = l.s. distance, 
#   but without dx and dy !!!! has to be changed if ever needed
# @details raw version that corresponds to spatstat.closepairs in call
#    but not in safety.
#    if invscale is given, it is used for the calculation, otherwise lambda is used
#    if none of them is given, the euclidean distances are returned
#' @rdname sosspp-internal
#' @keywords internal
#' @export
#' 
lsclosepairs <- function (X, rmax, invscale, lambda, 
               ordered = TRUE, what = c("all", "indices"))
{
  use.scale <- !missing(invscale)
  use.lambda <- !use.scale && !missing(lambda)
  use.none <- !use.scale & !use.lambda
  what <- match.arg(what)
  
  if (use.scale)  dists <- lspairdist(X$x, X$y, invscale = invscale)
  else if (use.lambda) dists <- lspairdist(X$x, X$y, lambda = lambda)
  else dists <- pairdist(X$x, X$y)
  
  if (use.none) {
    lsdist <- dists
    eudist <- dists
  } else {
    lsdist <- dists$lsd
    eudist <- dists$eud
  }  
  indi <- which(lsdist <= rmax, arr.ind = TRUE)
  i <- indi[, "row"]
  j <- indi[, "col"]
  OK <- i != j
  i <- i[OK]
  j <- j[OK]
  answer <- switch(what, 
            all = list(i = i, j = j, 
                  xi = X$x[i], yi = X$y[i], xj = X$x[j],  yj = X$y[j],
                  # dx = dx, dy = dy, 
                  d = lsdist[indi[OK, ]],
                  eud = eudist[indi[OK, ]]),
           indices = list(i = i, j = j))
  return(answer)
}