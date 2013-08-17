#' @useDynLib hidden2ndorder
# @param x,y vectors of x and y coordinates
# @param la intensity vector
#' @rdname hidden2ndorder-internal
#' @keywords internal
#' @export
#'


lspairdist <- function(x, y, la)
{
   .Call("lspairdist", x, y, la)
}
  

# @param X point pattern
# @param lambda a vector of intensities in the points of X
# @param rmax maximum locally scaled distance 
# @param ordered logical, is ignored
# @returns the locally scaled close pairs, where d= Euclidean distance, 
#   lsd = l.s. distance, 
#   but without dx and dy !!!! has to be changed if ever needed
# @details raw version that corresponds to spatstat.closepairs in call
#    but not in safety.
#' @rdname hidden2ndorder-internal
#' @keywords internal
#' @export
#' 
lsclosepairs <- function (X, lambda, rmax, ordered = TRUE, what = c("all", "indices"))
{
  what <- match.arg(what)
  dists <- lspairdist(X$x, X$y, lambda)
  lsdist <- dists$lsd
  eudist <- dists$eud
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