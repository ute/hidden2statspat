########################################################################
#                                                                      #
# Analysis of the Scholtzia Data Set as in Hahn & Jensen (2013)        #
#                                                                      #
########################################################################

# load the scholtzia data set and the intensity estimates used in H&J(2013) from sostatpp

data(scholtzia)
data(intensities)

# the data show higher intensity at the bottom than at the top of the plot,
# therefore we divide it into two halves according to y-coordinate

mediany <- median(scholtzia$y)
rangey <- scholtzia$window$yrange
rangex <- scholtzia$window$yrange

testset <- tiles(quadrats(bronzefilter, xbreaks = rangex, 
                 ybreaks = c(rangey[1], mediany, rangey[2])))

quads1 <- tiles(quadrats(testset[[2]], nx = 4, ny = 1))
quads2 <- tiles(quadrats(testset[[1]], nx = 2, ny = 2))

# colors for plotting, and a plot of the quadrats (Figure 19)

style1 <- list(col = "red", light = .6)
style2 <- list(col = "blue", light = .6)

quadratsplot(scholtzia, quads1, quads2, style1, style2, pch=16, cex=.5)

##### ----------  analysis as locally rescaled second-order stationary  --------

scholtzia.s <- rescaled(scholtzia, lambda = scholtzia.intens)

Kpermute.test(scholtzia.s, quads1 = quads1, quads2 = quads2, rmax = 1.25, 
              use.tbar = TRUE)

# visualisation: Figure 20, upper left

plot(estOnQuadrats(scholtzia.s, fun = estK, quads = quads1, rmax = 1.25), 
     style = style1, ylim = c(0, 6))
plot(estOnQuadrats(scholtzia.s, fun = estK, quads = quads2, rmax = 1.25), 
     style = style2, add = TRUE)

##### ----------- adjusting intensity on quadrats by "normpower" ----------
Kpermute.test(scholtzia.s, quads1 = quads1, quads2 = quads2, rmax = 1.25, 
              use.tbar = TRUE, normpower = 2)

# visualisation: Figure 20, upper right

plot(estOnQuadrats(scholtzia.s, fun = estK, quads = quads1, rmax = 1.25,
                  normpower = 2), 
     style = style1, ylim = c(0, 6))
plot(estOnQuadrats(scholtzia.s, fun = estK, quads = quads2, rmax = 1.25,
                  normpower = 2), 
     style = style2, add = TRUE)


##### ----------  analysis as reweighted second-order stationary  --------

scholtzia.w <- reweighted(scholtzia, lambda = scholtzia.intens)

Kpermute.test(scholtzia.w, quads1 = quads1, quads2 = quads2, rmax = 2.0, 
              use.tbar = TRUE)

# visualisation: Figure 20, lower left

plot(estOnQuadrats(scholtzia.w, fun = estK, quads = quads1, rmax = 2.0), 
     style = style1, ylim = c(0, 15))
plot(estOnQuadrats(scholtzia.w, fun = estK, quads = quads2, rmax = 2.0), 
     style = style2, add = TRUE)


##### ----------- adjusting intensity on quadrats by "normpower" ----------


Kpermute.test(scholtzia.w, quads1 = quads1, quads2 = quads2, rmax = 2.0, 
              use.tbar = TRUE, normpower = 2)

# visualisation: Figure 20, lower right

plot(estOnQuadrats(scholtzia.w, fun = estK, quads = quads1, rmax = 2.0, normpower = 2), 
     style = style1, ylim = c(0, 15))
plot(estOnQuadrats(scholtzia.w, fun = estK, quads = quads2, rmax = 2.0, normpower = 2), 
     style = style2, add = TRUE)



