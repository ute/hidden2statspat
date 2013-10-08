########################################################################
#                                                                      #
# Analysis of the Bronze Filter Data Set as in Hahn & Jensen (2013)    #
#                                                                      #
########################################################################

data(bronzefilter)

# intensity estimate assuming a gradient in x direction, 
# with spatstat function rhohat

lambdarho <- rhohat(unmark(bronzefilter), function(x, y) x, bw = 1)
lambda <- predict(lambdarho)

#####----------  analysis as rescaled or reweighted ---------- 

# determine two sets of 3x3 quadrats, with approximately the same number
# of points, using quantiles of the x-values 

xmedian <- median(bronzefilter$x)
xbreaks <- quantile(c(bronzefilter$x, bronzefilter$window$xrange), (0:6)/6)
quads <- tiles(quadrats(bronzefilter, xbreaks = xbreaks, ny = 3))

group1 <- sapply(quads, function(w) min(w$x <= xmedian))
quads1 <- quads[group1 == 1]
quads2 <- quads[group1 == 0]

# Figure 10: plot of the test sets

# colors for plotting
style1 <- list(col = "red", light = .6)
style2 <- list(col = "blue", light = .6)

# could also be set to black and white printing, e.g.
# style1 <- list(col = "black", lty = "solid", light = .6)
# style2 <- list(col = "black", lty = "dashed", light = .8)

quadratsplot(bronzefilter, quads1, quads2, style1, style2, pch=16, cex=.5)


##### ----------  analysis as rescaled second-order stationary  --------

bronzs <- rescaled(bronzefilter, lambda = lambda)

# testing the hypothesis that the pattern is rescaled second-order stationary
#
# !!!! to save time, we set the number of permutations to 1000
# This leads to random test resuls
# let nperm = 25000 in order to obtain exactly the same results as in H&J
#

nperm <- 1000

Kpermute.test(bronzs, quads1 = quads1, quads2 = quads2, rmax = 1.25, 
              nperm = nperm)

# using the test statistic $\bar T$ instead

Kpermute.test(bronzs, quads1 = quads1, quads2 = quads2, rmax = 1.25, 
              use.tbar = TRUE, nperm = nperm)

# visualisation: Figure 11, left

plot(estOnQuadrats(bronzs, fun = estK, quads = quads1, rmax = 1.25), 
     style = style1, ylim = c(0, 5))
plot(estOnQuadrats(bronzs, fun = estK, quads = quads2, rmax = 1.25), 
     style = style2, add = TRUE)

##### ----------  analysis as reweighted second-order stationary  --------

bronzw <- reweighted(bronzefilter, lambda = lambda)

Kpermute.test(bronzw, quads1 = quads1, quads2 = quads2, rmax = 0.4,
              nperm = nperm)

# using the test statistic $\bar T$ instead

Kpermute.test(bronzw, quads1 = quads1, quads2 = quads2, rmax = 0.4, 
              use.tbar = TRUE, nperm = nperm)

# visualisation: Figure 11, right

plot(estOnQuadrats(bronzw, fun = estK, quads = quads1, rmax = 0.4), 
     style = style1, ylim = c(0, 0.5))
plot(estOnQuadrats(bronzw, fun = estK, quads = quads2, rmax = 0.4), 
     style = style2, add = TRUE)


##### ----------  analysis as retransformed second-order stationary  --------

# we use the fact that the transformation follows a gradient in x - direction,
# and apply the estimated intensity

bronzt <- retransformed(bronzefilter, trafo="gradx", lambda = lambda)

# as in the paper, we specify the subwindows on the backtransformed pattern
# note that this only necessary for obtaining a plot as Figure 12, but the
# test could also be carried out using the original quadrats as in the analysis
# of rescaled s.o.s.

# two testsets of same size
testset <- tiles(quadrats(bronzefilter, nx = 2, ny = 1))

quadt1 <- quadrats(testset[[1]], 3, 3)
quadt2 <- quadrats(testset[[2]], 3, 3)

# Figure 12, left part

quadratsplot(backtransformed(bronzt), quadt1, quadt2, style1, style2, 
             pch=16, cex=.5)

# Figure 12, right part

plot(estOnQuadrats(bronzt, fun = estK, tquads = quadt1, rmax = 0.5), 
     style = style1, ylim = c(0, 0.8))
plot(estOnQuadrats(bronzt, fun = estK, tquads = quadt2, rmax = 0.5), 
     style = style2, add = TRUE)

# testing

Kpermute.test(bronzt, tquads1 = quadt1, tquads2 = quadt2, rmax=.5, 
              nperm = nperm)
Kpermute.test(bronzt, tquads1 = quadt1, tquads2 = quadt2, rmax=.5, 
              nperm = nperm, use.tbar = TRUE)


#### directional analysis of the backtransformed pattern ---------------------

# Figure 14, left part, and test

plot(estOnQuadrats(bronzt, fun = estDeltaKdir, tquads = quadt1, rmax = 0.6), 
     style = style1, ylim = c(-.5, .5))
plot(estOnQuadrats(bronzt, fun = estDeltaKdir, tquads = quadt2, rmax = 0.6), 
     style = style2, add = TRUE)

Kpermute.test(bronzt, Kfun = estDeltaKdir, tquads1 = quadt1, tquads2 = quadt2, 
              rmax=.5, use.tbar = FALSE, nperm = nperm)

Kpermute.test(bronzt, Kfun = estDeltaKdir, tquads1 = quadt1, tquads2 = quadt2, 
              rmax=.5, use.tbar = TRUE, nperm = nperm)


#### directional analysis of the locally rescaled pattern ---------------------

# Figure 14, right part, and test

plot(estOnQuadrats(bronzs, fun = estDeltaKdir, quads = quads1, rmax = 1.25), 
     style = style1, ylim = c(-3, 3))
plot(estOnQuadrats(bronzs, fun = estDeltaKdir, quads = quads2, rmax = 1.25), 
     style = style2, add = TRUE)

Kpermute.test(bronzs, Kfun = estDeltaKdir, quads1 = quads1, quads2 = quads2, 
              rmax=1.25, use.tbar = FALSE, nperm = nperm)

Kpermute.test(bronzs, Kfun = estDeltaKdir, quads1 = quads1, quads2 = quads2, 
              rmax=1.25, use.tbar = TRUE, nperm = nperm)


