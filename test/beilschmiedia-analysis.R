# analysis of beilschmiedia

data(bei)
data(intensities)

style1 <- list(col = "red",  light = .6)
style2 <- list(col = "blue", light = .6)

allquads <- tiles(quadrats(bei, nx=8, ny=4))
npts <- sapply(allquads, function(w) npoints(bei[w])) # replaced quadratcount

enoughpoints <- npts > 30

medin <- median(npts[enoughpoints])
few <- npts <= medin & enoughpoints
many <- npts >= medin

#allquads <- tiles(allquads)

quads1 <- allquads[many]
quads2 <- allquads[few]
quads0 <- allquads[!enoughpoints]

quadratsplot(bei, quads1, quads2, style1, style2, quads0 = quads0, pch=16, cex=.3)



Kpermute.test(bei, type = "w", lambda=bei.intens.bayes, 
              quads1 = quads1, quads2 = quads2, rmax = 25, nperm = 1000)

#   Studentized permutation test of reWeighted hidden second-order stationarity,
# 	isotropic version, using K_0
# 	test statistic: T, integration limits 0 to 25
# 	exact test, using all permutations
# 
# data:  point pattern bei
# T = 0.0366, p-value = 0.001113
# alternative hypothesis: not reWeighted second-order stationary 

Kpermute.test(bei, type = "w", lambda=bei.intens.bayes, 
              quads1 = quads1, quads2 = quads2, rmax = 25, use.tbar=TRUE, nperm = 1000)

#   Studentized permutation test of reWeighted hidden second-order stationarity,
# 	isotropic version, using K_0
# 	test statistic: Tbar, integration limits 0 to 25
# 	using 25000 randomly selected permuations
# 
# data:  point pattern bei
# Tbar = 0.038, p-value = 0.00116
# alternative hypothesis: not reWeighted second-order stationary 

Kpermute.test(bei, type = "w", lambda=bei.intens.maxlik, 
              quads1 = quads1, quads2 = quads2, rmax = 25, use.tbar=TRUE, nperm = 1000)

#   Studentized permutation test of reWeighted hidden second-order stationarity,
# 	isotropic version, using K_0
# 	test statistic: Tbar, integration limits 0 to 25
# 	using 25000 randomly selected permuations
# 
# data:  point pattern bei
# Tbar = 0.0297, p-value = 0.00024
# alternative hypothesis: not reWeighted second-order stationary 

Kw1 <- estOnQuadrats(bei, type="w", lambda=bei.intens.bayes, quads = quads1, rmax = 25)
Kw2 <- estOnQuadrats(bei, type="w", lambda=bei.intens.bayes, quads = quads2, rmax = 25)
plot(Kw1, col = 1, lty = "solid", ylim=c(0,3000))
plot(Kw2, col = 1, lty = "dashed", add=T)

Kw1 <- estOnQuadrats(bei, type="w", lambda=bei.intens.maxlik, quads = quads1, rmax = 25)
Kw2 <- estOnQuadrats(bei, type="w", lambda=bei.intens.maxlik, quads = quads2, rmax = 25)
plot(Kw1, col = 1, lty = "solid", ylim=c(0,3000))
plot(Kw2, col = 1, lty = "dashed", add=T)

Kpermute.test(bei, type = "w", lambda=bei.intens.nonpar, 
              quads1 = quads1, quads2 = quads2, rmax = 25, use.tbar=TRUE, nperm = 1000)

#   Studentized permutation test of reWeighted hidden second-order stationarity,
# 	isotropic version, using K_0
# 	test statistic: Tbar, integration limits 0 to 25
# 	using 25000 randomly selected permuations
# 
# data:  point pattern bei
# Tbar = 0.0099, p-value = 0.1632
# alternative hypothesis: not reWeighted second-order stationary 

Kpermute.test(bei, type = "w", lambda=bei.intens.maxlik, normpower = 2, 
              quads1 = quads1, quads2 = quads2, rmax = 25, use.tbar=TRUE, nperm=1000)
#   Studentized permutation test of reWeighted hidden second-order stationarity,
# 	isotropic version, using K_0
# 	test statistic: Tbar, integration limits 0 to 25
# 	using 5000 randomly selected permuations
# 
# data:  point pattern bei
# Tbar = 6e-04, p-value = 0.8748
# alternative hypothesis: not reWeighted second-order stationary 

Kpermute.test(bei, type = "w", lambda=bei.intens.bayes, normpower = 2, 
              quads1 = quads1, quads2 = quads2, rmax = 25, use.tbar=TRUE, nperm=1000)

#   Studentized permutation test of reWeighted hidden second-order stationarity,
# 	isotropic version, using K_0
# 	test statistic: Tbar, integration limits 0 to 25
# 	using 5000 randomly selected permuations
# 
# data:  point pattern bei
# Tbar = 0.001, p-value = 0.8684
# alternative hypothesis: not reWeighted second-order stationary 
