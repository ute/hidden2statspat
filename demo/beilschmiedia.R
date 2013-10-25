########################################################################
#                                                                      #
# Analysis of the Beilschmiedia Data Set as in Hahn & Jensen (2013)    #
#                                                                      #
########################################################################

# load the beilschmiedia data from spatstat, 
# and the intensity estimates used in H&J(2013) from sostatpp

data(bei)
data(intensities)

# partition the point pattern into 8 x 4 quadrats, and select subsets 
# with high and low intensity

beiquads <- twoquadsets(bei, nx = 8, ny = 4, minpoints = 30)

# generate a plot of the quadrats
# plot style for the two used subsets

beistyle <- list(hi = simplist(col="red",  alpha=.4, col.win="red", alpha.win=.4),
                 lo = simplist(col="blue", alpha=.4, col.win="blue"))

quadratsplot(bei, beiquads, beistyle, pch = 16, cex = .4, 
  main = "Beilschmiedia quadrats for testing")

# testing reweighted second-order stationarity ----------------------------

# The number of permutations is here set to 1000 only, to save time. 
# In the paper, we used nperm = 5000. Since the tests are random, the results
# may vary slightly

nperm <- 1000

# tests based on test statistic \bar{T} 

# tagging the beilschmiedia data set as reweighted, with Bayesian intensity estimate

bei.wbay <- reweighted(bei, intensity = bei.intens.bayes)

testresult.wbay <- sos.test(bei.wbay, beiquads, rmax = 25, use.tbar=TRUE, nperm = nperm)

testresult.wbay

#plot Figure 17, left

plot(testresult.wbay, beistyle, ylim = c(0,3000))


# using the maximum likelihood intensity estimate

# tagging the beilschmiedia data set as reweighted, with Bayesian intensity estimate

bei.wml <- reweighted(bei, intensity = bei.intens.maxlik)

testresult.wml <- sos.test(bei.wml, beiquads, rmax = 25, use.tbar=TRUE, nperm = nperm)

testresult.wml

#plot Figure 17, right

plot(testresult.wml, beistyle, ylim = c(0,3000))

####### same, with nonparametric intensity estimate, 

sos.test(reweighted(bei, intensity = bei.intens.nonpar),
              beiquads, rmax = 25, 
              use.tbar=TRUE, nperm = nperm)

####### same, with parameter normpower = 2

sos.test(bei.wml, beiquads, rmax = 25, normpower = 2,
              use.tbar=TRUE, nperm = nperm)

sos.test(bei.wbay, beiquads, rmax = 25, normpower = 2,
              use.tbar=TRUE, nperm = nperm)

sos.test(reweighted(bei, intensity = bei.intens.nonpar),
              beiquads, rmax = 25, normpower = 2,
              use.tbar=TRUE, nperm = nperm)
