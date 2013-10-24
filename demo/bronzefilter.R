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
# of points, using quantiles of the x-values. Assign to the sets according to
# x -ccordinate

bronzequads <- twoquadsets(bronzefilter, nx = 6, ny = 3, gradx = TRUE, 
                           by.intensity = FALSE)

# Figure 10: plot of the test sets

# colors for plotting
styles <- list( lo = style(col="red", col.win="red", alpha.win=.4, alpha=.4), 
                hi = style(col="blue", col.win="blue", alpha.win=.4, alpha=.4))


# could also be set to black and white printing, e.g.
# style1 <- list(col = "black", lty = "solid", light = .6)
# style2 <- list(col = "black", lty = "dashed", light = .8)

quadratsplot(bronzefilter, bronzequads, styles, use.marks=FALSE, pch=16, cex=.5)


##### ----------  analysis as rescaled second-order stationary  --------

bronzs <- rescaled(bronzefilter, intensity = lambda)

# testing the hypothesis that the pattern is rescaled second-order stationary
#
# !!!! to save time, we set the number of permutations to 1000
# This leads to random test resuls
# let nperm = 25000 in order to obtain exactly the same results as in H&J
#

nperm <- 1000

test_s <- sos.test(bronzs, bronzequads, rmax = 1.25, nperm = nperm)
test_s

plot(test_s, styles)

# using the test statistic $\bar T$ instead

sos.test(bronzs, bronzequads, rmax = 1.25, nperm = nperm, use.tbar = TRUE)

# visualisation: Figure 11, left

plot(test_s, styles)

##### ----------  analysis as reweighted second-order stationary  --------

bronzw <- reweighted(bronzefilter, intensity = lambda)

test_w<- sos.test(bronzw, bronzequads, rmax = 0.4, nperm = nperm)
test_w

# using the test statistic $\bar T$ instead
sos.test(bronzw, bronzequads, rmax = 0.4, nperm = nperm, use.tbar = TRUE)

# visualisation: Figure 11, right

plot(test_w, styles)

##### ----------  analysis as retransformed second-order stationary  --------

# we use the fact that the transformation follows a gradient in x - direction,
# and apply the estimated intensity

bronzt <- retransformed(bronzefilter, backtrafo="gradx", intensity = lambda)

test_t<- sos.test(bronzt, bronzequads, rmax = 0.5, nperm = nperm)
test_t

# Figure 12, left part: backtransformed quadratsplot 
quadratsplot(bronzt, bronzequads, styles, use.marks = FALSE, pch=16, cex=.5,
  backtransformed = TRUE)

# Figure 12, right part

plot(test_t, styles)

#### directional analysis of the backtransformed pattern ---------------------

test_tdir<- sos.test(bronzt, bronzequads, Kfun = DeltaKdir.est, rmax = 0.5, nperm = nperm)
print(test_tdir)

test_sdir<- sos.test(bronzs, bronzequads, Kfun = DeltaKdir.est, rmax = 1.25, nperm = nperm)
print(test_sdir)

# Figure 14, left part, and test

plot(test_tdir, styles)

plot(test_sdir, styles)
