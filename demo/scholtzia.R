########################################################################
#                                                                      #
# Analysis of the Scholtzia Data Set as in Hahn & Jensen (2013)        #
#                                                                      #
########################################################################

# load the scholtzia data set and the intensity estimates used in H&J(2013) from sostatpp

data(scholtzia)
data(intensities)

# "home made" quadrats:
#
# the data show higher intensity at the bottom than at the top of the plot,
# therefore we divide it into two halves according to y-coordinate

testset <- twoquadsets(scholtzia, nx = 1, ny = 2, grady = TRUE)

# now the subsets are subdivided into 2x2 and 4x1 quadrats, respectively
# new lo: low intensity, the upper part (was hi), 2x2  hi: lower part, 4x1

quadsets <- list(lo = quadrats(testset$hi[[1]], nx = 2, ny = 2),
                 hi = quadrats(testset$lo[[1]], nx = 4, ny = 1))

# colors for plotting, and a plot of the quadrats (Figure 19)

styles <-  list(hi = simplist(col="red",  alpha=.4, col.win="red", alpha.win=.4),
               lo = simplist(col="blue", alpha=.4, col.win="blue"))

quadratsplot(scholtzia, quadsets, styles, pch=16, cex=.5)

##### ----------  analysis as locally rescaled second-order stationary  --------

scholtzia_s <- rescaled(scholtzia, intensity = scholtzia.intens)

test_s <- sos.test(scholtzia_s, quadsets, rmax = 1.25, use.tbar = TRUE)
print(test_s)
# visualisation: Figure 20, upper left
plot(test_s, styles)

##### ----------- adjusting intensity on quadrats by "normpower" ----------

test_s2 <- sos.test(scholtzia_s, quadsets, rmax = 1.25, use.tbar = TRUE, normpowe = 2)
print(test_s2)
# visualisation: Figure 20, upper right
plot(test_s2, styles)


##### ----------  analysis as locally rescaled second-order stationary  --------

scholtzia_w <- reweighted(scholtzia, intensity = scholtzia.intens)

test_w <- sos.test(scholtzia_w, quadsets, rmax = 2, use.tbar = TRUE)
print(test_w)
# visualisation: Figure 20, lower left
plot(test_w, styles)

##### ----------- adjusting intensity on quadrats by "normpower" ----------

test_w2 <- sos.test(scholtzia_w, quadsets, rmax = 2, use.tbar = TRUE, normpower = 2)
print(test_w2)
# visualisation: Figure 20, lower right
plot(test_w2, styles)

