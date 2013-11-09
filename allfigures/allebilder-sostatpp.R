# all images for our paper

library(sostatpp)

ploptions(alphamix = TRUE)


setparams <- function(bw = FALSE) {
  if (bw)
  {  
    styles <<- list( 
      hi = simplist(col = "black", lwd = 0.8, lty = "solid", col.win = gray(.2), 
                     alpha = .5, alpha.win = .3),
      lo = simplist(col = gray(.2), lwd = 0.8, lty = "dashed", col.win = gray(.5), 
                     alpha = .5, alpha.win = .3)
      )
    theostyle <<- simplist(lty = "dotted", lwd = 1.5)
  } else {
    col.rot  <<- rgb(1, 0, .1) #"red",
    col.blau <<- rgb(0, .2, 1) #"blue",
    styles <<- list( 
      hi = simplist(col = col.rot, lwd = 0.8, lty = "solid", col.win = col.rot, 
                     alpha = .4, alpha.win = .3, alpha.smf = .8),
      lo = simplist(col = col.blau, lwd = 0.8, lty = "solid", col.win = col.blau, 
                     alpha = .4, alpha.win = .3, alpha.smf = .8)
      )
    theostyle <<- simplist(lty = "dashed", lwd = 1.5)
  }
}

patterns <- function() par(mar=c(0.01, 0.01, 0.01, 0.01))
curves <- function() par(mar=c(3.6, 4.1, 1.1, 1.1)) 

setparams()

# bronze: hi and lo refer to x-coordinate. Therefore, swap colors
brostyles <- list(hi = styles$lo, lo = styles$hi)

data(intensities)


# BRONZE ------------------------------------------------------------------

bronzequads <- twoquadsets(bronzefilter, nx = 6, ny = 3, gradx = TRUE, 
                           by.intensity = FALSE)

pdf(file="bronz-quads.pdf", width=7, height = 3)
patterns()
quadratsplot(bronzefilter, bronzequads, brostyles, use.marks=FALSE, pch=16, cex=.5,
  main = "")
dev.off()


bronze.w <- reweighted(bronzefilter, intensity = bronze.intens)
bronze.s <- rescaled(bronzefilter, intensity = bronze.intens)
test_s <- sos.test(bronze.s, bronzequads, rmax = 1.25, nperm = 2)
test_w<- sos.test(bronze.w, bronzequads, rmax = 0.4, nperm = 2)

pdf(file="bronze-Ks.pdf", width = 4, height = 3.2)
 curves()
 plot(test_s, brostyles, theostyle, ylim = c(0, 5))
dev.off()

pdf(file="bronze-Kw.pdf", width = 4, height = 3.2)
 curves()
 plot(test_w, brostyles, theostyle, ylim = c(0, .5))
dev.off()

bronze.t <- retransformed(bronzefilter, backtrafo="gradx", 
                          intensity = bronze.intens)

test_t<- sos.test(bronze.t, bronzequads, rmax = 0.5, nperm = 2)

# Figure 12, left part: backtransformed quadratsplot 
pdf(file="bronztra-quads.pdf", width=7, height = 3)
patterns()
quadratsplot(bronze.t, bronzequads, brostyles, use.marks = FALSE, pch=16, cex=.5,
  backtransformed = TRUE, main = "")
dev.off()


pdf(file="bronze-Kt.pdf", width = 4, height = 3.2)
curves()
 plot(test_t, brostyles, theostyle, ylim = c(0, .8))
dev.off()

test_tdir<- sos.test(bronze.t, bronzequads, Kfun = DeltaKdir.est, rmax = 0.5, nperm = 2)

test_sdir<- sos.test(bronze.s, bronzequads, Kfun = DeltaKdir.est, rmax = 1.25, nperm = 2)

pdf(file="bronze-Kdiff.pdf", width = 4, height = 3.2)
 curves()
 plot(test_sdir, brostyles, theostyle, ylim = c(-3, 3))
dev.off()

pdf(file="bronztra-Kdiff.pdf", width = 4, height = 3.2)
 curves()
 plot(test_tdir, brostyles, theostyle, ylim = c(-.5, .5))
dev.off()


# BEILSCHMIEDIA -----------------------------------------------------------

data(bei)
data(intensities)

q8x4 <- quadrats(bei, nx=8, ny=4)


nco <- 4096
#colseq <- gray(seq(0.00001,1,len=4096)^.2)
colseq <- gray(seq(0.000001,1,len=nco))
colseq <- colseq[length(colseq):1]
refr <- log(range(bei.intens.nonpar$v))
minr <- min(refr)
rbay <- round((nco-1)*((log(range(bei.intens.bayes$v))-minr)/diff(refr))+1)
rml <- round((nco-1)*((log(range(bei.intens.maxlik$v))-minr)/diff(refr))+1)

beiqs <- function() par(mar=c(0.1, 0.1, 0.1, 2.1))

pdf("bei-qintens-bayes.pdf", width=6, height = 2.5)
beiqs()
plot(bei.intens.bayes, log=T, main="", col=colseq[rbay[1]:rbay[2]])
plot(q8x4, add=T)
dev.off()

pdf("bei-qintens-maxlik.pdf", width=6, height = 2.5)
beiqs()
plot(bei.intens.maxlik, log=T, main="", col=colseq[rml[1]:rml[2]])
plot(q8x4, add=T)
dev.off()

pdf("bei-qintens-nonpar.pdf", width=6, height = 2.5)
beiqs()
plot(bei.intens.nonpar, log=T, main="", col=colseq)
plot(q8x4, add=T)
dev.off()


beiquads <- twoquadsets(bei, nx = 8, ny = 4, minpoints = 30)

pdf(file="bei-quads-npoints-4.pdf", width=8, height = 4)
patterns()
quadratsplot(bei, beiquads, styles, pch = 16, cex = .4, main = "")
dev.off()


bei.wbay <- reweighted(bei, intensity = bei.intens.bayes)
bei.wml <- reweighted(bei, intensity = bei.intens.maxlik)
bei.wnon <- reweighted(bei, intensity = bei.intens.nonpar)

test.wbay <- sos.test(bei.wbay, beiquads, rmax = 25, use.tbar=TRUE, nperm = 2)
test.wml <- sos.test(bei.wml, beiquads, rmax = 25, use.tbar=TRUE, nperm = 2)
test.wnon <- sos.test(bei.wnon, beiquads, rmax = 25, use.tbar=TRUE, nperm = 2)
test.wbay2 <- sos.test(bei.wbay, beiquads, rmax = 25, use.tbar=TRUE, normpower = 2,nperm = 2)
test.wml2 <- sos.test(bei.wml, beiquads, rmax = 25, use.tbar=TRUE, normpower = 2,nperm = 2)
test.wnon2 <- sos.test(bei.wnon, beiquads, rmax = 25, use.tbar=TRUE, normpower = 2, nperm = 2)


#plot Figure 17, left
pdf(file="bei-Kw-bayes.pdf", width = 4, height = 4)
curves()
 plot(test.wbay, styles, theostyle, ylim = c(0, 3000))
dev.off()

pdf(file="bei-Kw-maxli.pdf", width = 4, height = 4)
curves()
 plot(test.wml, styles, theostyle, ylim = c(0, 3000))
dev.off()

pdf(file="bei-Kw-nonpar.pdf", width = 4, height = 4)
curves()
 plot(test.wnon, styles, theostyle, ylim = c(0, 3000))
dev.off()

pdf(file="bei-Kw-bayes-adj.pdf", width = 4, height = 4)
curves()
 plot(test.wbay2, styles, theostyle, ylim = c(0, 3000))
dev.off()

pdf(file="bei-Kw-maxli-adj.pdf", width = 4, height = 4)
curves()
 plot(test.wml2, styles, theostyle, ylim = c(0, 3000))
dev.off()

pdf(file="bei-Kw-nonpar-adj.pdf", width = 4, height = 4)
curves()
 plot(test.wnon2, styles, theostyle, ylim = c(0, 3000))
dev.off()


# SCHOLTZIA median --------------------------------------------------------

data(scholtzia)

# schiefe quadrates
W1 <- owin(c(0, 22), c(0, 4.4))
W2 <-  owin(c(0, 22), c(4.4, 22))

quads1 <- tiles(quadrats(W1, nx=4, ny=1))
quads2 <- tiles(quadrats(W2, nx=2, ny=2))

quads1[[1]] <- owin(poly=list(x=c(0, 5.5, 5.5, 0), y=c(0, 0, 5.2, 6.0)))
quads1[[2]] <- owin(poly=list(x=c(5.5,11, 11, 5.5), y=c(0, 0, 4.4, 5.2)))
quads2[[3]] <- owin(poly=list(x=c(0, 11, 11, 0), y=c(6.0, 4.4, 13.2, 13.2)))

# ordentliche quadrats
testset <- twoquadsets(scholtzia, nx = 1, ny = 2, grady = TRUE)

# now the subsets are subdivided into 2x2 and 4x1 quadrats, respectively
# new lo: low intensity, the upper part (was hi), 2x2  hi: lower part, 4x1

quadsets <- list(lo = quadrats(testset$hi[[1]], nx = 2, ny = 2),
                 hi = quadrats(testset$lo[[1]], nx = 4, ny = 1))

# colors for plotting, and a plot of the quadrats (Figure 19)

pdf(file="scholtzia-quads.pdf", width=4, height = 4)
patterns()
quadratsplot(scholtzia, quadsets, styles, pch = 16, cex=.5, main = "")
dev.off()

scholtzia_s <- rescaled(scholtzia, intensity = scholtzia.intens)
scholtzia_w <- reweighted(scholtzia, intensity = scholtzia.intens)

test_s <- sos.test(scholtzia_s, quadsets, rmax = 1.25, use.tbar = TRUE)
test_s2 <- sos.test(scholtzia_s, quadsets, rmax = 1.25, normpower = 2, use.tbar = TRUE)
test_w <- sos.test(scholtzia_w, quadsets, rmax = 2, use.tbar = TRUE)
test_w2 <- sos.test(scholtzia_w, quadsets, rmax = 2, normpower = 2, use.tbar = TRUE)

pdf(file="scholtzia-Kw-orig.pdf", width = 4, height = 3.2)
curves()
 plot(test_w, styles, theostyle, ylim = c(0, 15))
dev.off()

pdf(file="scholtzia-Kw-adj.pdf", width = 4, height = 3.2)
curves()
 plot(test_w2, styles, theostyle, ylim = c(0, 15))
dev.off()

pdf(file="scholtzia-Ks-orig.pdf", width = 4, height = 3.2)
curves()
 plot(test_s, styles, theostyle, ylim = c(0, 6))
dev.off()

pdf(file="scholtzia-Ks-adj.pdf", width = 4, height = 3.2)
curves()
 plot(test_s2, styles, theostyle, ylim = c(0, 6))
dev.off()