
#####################################################################
#
#        edgecrossRipley.R
#
#    Ripley isotropic edge correction weights, cross directional 
#
#    derived from spatstat's Ripley edge correction
#
#  edgecross.Ripley(X, r, dphi, W)      compute isotropic correction weights
#                            for centres X[i], radii r[i,j], window W, 
#                  directions 0+- dphi, pi/2+- dphi
#
#
#######################################################################

# @param X point pattern of class ppp
# @param r matrix of radii
# @param W window
# @param dphi half angle
#' @rdname sostatpp-internal
#' @keywords internal
#' @export

edgecross.Ripley <- local({
  
  small <- function(x) { abs(x) < .Machine$double.eps }
  
  hang <- function(d, r) {
    nr <- nrow(r)
    nc <- ncol(r)
    answer <- matrix(0, nrow=nr, ncol=nc)
    # replicate d[i] over j index
    d <- matrix(d, nrow=nr, ncol=nc)
    hit <- (d < r)
    answer[hit] <- acos(d[hit]/r[hit])
    answer
  }
  
  edgecross.Ripley <- function(X, r, dphi = pi/4, W=X$window,  maxweight=100) {
    # X is a point pattern, or equivalent
    X <- as.ppp(X, W)
    W <- X$window
    
    switch(W$type,
           rectangle={},
           polygonal={
              stop("sorry, Ripley isotropic cross correction only implemented for rectangular windows")
           },
           mask={
             stop("sorry, Ripley isotropic cross correction only implemented for rectangular windows")
           }
    )
    
    n <- npoints(X)
    
    if(is.matrix(r) && nrow(r) != n)
      stop("the number of rows of r should match the number of points in X")
    if(!is.matrix(r)) {
      if(length(r) != n)
        stop("length of r is incompatible with the number of points in X")
      r <- matrix(r, nrow=n)
    }
    
    #
    Nr <- nrow(r)
    Nc <- ncol(r)
    if(Nr * Nc == 0) return(r)
    
    ##########
    
    x <- X$x
    y <- X$y
    
    ######## interpreted R code for rectangular case #########
    
    # perpendicular distance from point to each edge of rectangle
    # L = left, R = right, D = down, U = up
    dL  <- x - W$xrange[1]
    dR  <- W$xrange[2] - x
    dD  <- y - W$yrange[1]
    dU  <- W$yrange[2] - y
    
    # detect whether any points are corners of the rectangle
    corner <- (small(dL) + small(dR) + small(dD) + small(dU) >= 2)
    
    # angle between (a) perpendicular to edge of rectangle
    # and (b) line from point to corner of rectangle
    bLU <- atan2(dU, dL)
    bLD <- atan2(dD, dL)
    bRU <- atan2(dU, dR)
    bRD <- atan2(dD, dR)
    bUL <- atan2(dL, dU)
    bUR <- atan2(dR, dU)
    bDL <- atan2(dL, dD)
    bDR <- atan2(dR, dD)
    
    # The above are all vectors [i]
    # Now we compute matrices [i,j]
    
    # half the angle subtended by the intersection between
    # the circle of radius r[i,j] centred on point i
    # and each edge of the rectangle (prolonged to an infinite line)
    
    aL <- hang(dL, r)
    aR <- hang(dR, r)
    aD <- hang(dD, r) 
    aU <- hang(dU, r)
    
    # apply maxima
    # note: a* are matrices; b** are vectors;
    # b** are implicitly replicated over j index
    cLU <- pmin(aL, bLU)
    cUL <- pmin(aU, bUL)
    cLD <- pmin(aL, bLD)
    cDL <- pmin(aD, bDL)
    cRU <- pmin(aR, bRU)
    cUR <- pmin(aU, bUR)
    cRD <- pmin(aR, bRD)
    cDR <- pmin(aD, bDR)
    
    # interior angles
    intRU <- pmin(dphi, pi/2 - cUR) - pmin (dphi, cRU)
    intUR <- pmin(dphi, pi/2 - cRU) - pmin (dphi, cUR)
    intLU <- pmin(dphi, pi/2 - cUL) - pmin (dphi, cLU)
    intUL <- pmin(dphi, pi/2 - cLU) - pmin (dphi, cUL)
    intRD <- pmin(dphi, pi/2 - cDR) - pmin (dphi, cRD)
    intDR <- pmin(dphi, pi/2 - cRD) - pmin (dphi, cDR)
    intLD <- pmin(dphi, pi/2 - cDL) - pmin (dphi, cLD)
    intDL <- pmin(dphi, pi/2 - cLD) - pmin (dphi, cDL)
    
    intRL <- intRU + intRD + intLU + intLD
    intUD <- intUR + intUL + intDR + intDL
    
    # set conrers to dphi
    if(any(corner)){ intRL[corner,] <- dphi; intUD[corner,] <- dphi}
    
    # OK, now compute weight
    wRL <- 4*dphi / intRL
    wUD <- 4*dphi / intUD
    # eliminate wild values
    wUD <- matrix(pmax(0, pmin(maxweight, wUD)),
                     nrow=Nr, ncol=Nc)
    wRL <- matrix(pmax(0, pmin(maxweight, wRL)),
                  nrow=Nr, ncol=Nc)
    return(data.frame(w.horiz=wRL, w.vert=wUD))
  }    
  edgecross.Ripley
}) 