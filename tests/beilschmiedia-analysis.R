########################################################################
#                                                                      #
# Analysis of the Beilschmiedia Data Set as in Hahn & Jensen (2013)    #
#                                                                      #
########################################################################
#
# 07 october 2013, Ute Hahn (ute@imf.au.dk)
#
#
# The methods proposed in Hahn & Jensen (2013) are bundled in package
# sostatpp. It is available on R-forge, and depends on fdnonpar and spatstat.
# To run this code, first install spatstat from CRAN, then fnonpar and sostatpp
# from R-forge, by issuing the commands
#
# install.packages("spatstat")
# install.packages(c("fdnonpar", "sostatpp"), repos="http://R-Forge.R-project.org")
#
#

library(sostatpp)

# load the beilschmiedia data from spatstat, 
# and the intensity estimates used in H&J(2013) from sostatpp

data(bei)
data(intensities)

# partition the point pattern into quadrats, and select subsets 
# with high and low intensity

allquads <- tiles(quadrats(bei, nx=8, ny=4))
npts <- sapply(allquads, function(w) npoints(bei[w])) # replaced quadratcount

enoughpoints <- npts > 30

median_npts <- median(npts[enoughpoints])
few  <- npts <= median_npts & enoughpoints
many <- npts >= median_npts

# subsets of high and low intensity, and unused quadrats

quads1 <- allquads[many]
quads2 <- allquads[few]
quads0 <- allquads[!enoughpoints]

# plot style for the two used subsets

style1 <- list(col = "red",  light = .6)
style2 <- list(col = "blue", light = .6)

# a plot of the subsamples

quadratsplot(bei, quads1, quads2, style1, style2, 
             quads0 = quads0, pch=16, cex=.3)


# testing reweighted second-order stationarity ----------------------------

# The number of permutations is here set to 1000 only, to save time. 
# In the paper, we used nperm = 5000. Since the tests are random, the results
# may vary slightly

nperm <- 1000

# tests based on test statistic \bar{T} 

# tagging the beilschmiedia data set as reweighted, with Bayesian intensity estimate

bei.wbay <- reweighted(bei,  lambda=bei.intens.bayes)

Kpermute.test(bei.wbay, quads1 = quads1, quads2 = quads2, rmax = 25, 
              use.tbar=TRUE, nperm = nperm)

#plot Figure 17, left

plot(estOnQuadrats(bei.wbay, fun = estK, quads = quads1, rmax = 25), 
     style = style1, ylim = c(0, 3000))
plot(estOnQuadrats(bei.wbay, fun = estK, quads = quads2, rmax = 25), 
     style = style2, add = T)


# using the maximum likelihood intensity estimate

# tagging the beilschmiedia data set as reweighted, with Bayesian intensity estimate

bei.wml <- reweighted(bei, lambda=bei.intens.maxlik)

Kpermute.test(bei.wml, quads1 = quads1, quads2 = quads2, rmax = 25, 
              use.tbar=TRUE, nperm = nperm)

#plot Figure 17, right

plot(estOnQuadrats(bei.wml, fun = estK, quads = quads1, rmax = 25), 
     style = style1, ylim = c(0, 3000))
plot(estOnQuadrats(bei.wml, fun = estK, quads = quads2, rmax = 25), 
     style = style2, add = T)

####### same, with nonparametric intensity estimate, 

Kpermute.test(reweighted(bei, lambda = bei.intens.nonpar),
              quads1 = quads1, quads2 = quads2, rmax = 25, 
              use.tbar=TRUE, nperm = nperm)

####### same, with parameter normpower = 2

Kpermute.test(bei.wbay, quads1 = quads1, quads2 = quads2, rmax = 25, 
              normpower = 2,  
              use.tbar=TRUE, nperm = nperm)

Kpermute.test(bei.wml, quads1 = quads1, quads2 = quads2, rmax = 25, 
              normpower = 2,  
              use.tbar=TRUE, nperm = nperm)

