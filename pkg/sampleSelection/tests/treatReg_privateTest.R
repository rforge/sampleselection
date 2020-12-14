### Testing treatreg.
### These are tests that are not supposed to be included in CRAN

DGP <- function(N=1000, sigma=1, rho=0.8,
                alpha0=-1, alpha1=1, alpha2=1,
                beta0=0, beta1=1, beta2=1) {
   ## Generate random data
   library(mvtnorm)
   Sigma <- matrix(c(1, rho*sigma, rho*sigma, sigma^2), 2, 2)
   uv <- rmvnorm(N, mean=c(0,0), sigma=Sigma)
   u <- uv[,1]
   v <- uv[,2]
   x <- rnorm(N)
   z <- rnorm(N)
   ySX <- alpha0 + alpha1*x + alpha2*z + u
   yS <- ySX > 0
   yO <- beta0 + beta1*x + beta2*yS + v
   data.frame(yO, yS, x, z, ySX, u, v)
}

library(sampleSelection)
# the following command makes sure that sample() returns the same pseudo-random
# numbers in R 3.5.X and in R-devel  
suppressWarnings( RNGversion( "3.5.0" ) )
set.seed(1)
options(digits=3)
cat("NA, Inf in data.  Should show 93 observations\n")
dat <- DGP(100)
dat$yO[1] <- NA
dat$yO[2] <- Inf
dat$yS[3] <- NA
dat$x[5] <- NA
dat$x[6] <- Inf
dat$z[7] <- NA
dat$z[8] <- Inf
m <- treatReg(yS~x+z, yO~yS+x, data=dat)
print(summary(m))

mf <- model.frame(m)
print(mf[sample(nrow(mf), 10),])
mm <- model.matrix(m)
print(mm[sample(nrow(mm), 10),])
mm <- model.matrix(m, part="selection")
print(mm[sample(nrow(mm), 10),])

## Now test prediction
cat("predicted and actual selection values\n")
pl <- predict(m, part="selection", type="link")
pr <- predict(m, part="selection", type="response")
p <- cbind(pred.link=pl, "P[ys=1]"=pr, actual.resp=mf$yS)
print(p[sample(nrow(p), 10),])
cat("predicted and actual outcomes\n")
pu <- predict(m, part="outcome", type="unconditional")
pc <- predict(m, part="outcome", type="conditional")
p <- cbind("E[yo]"=pu, pc, "yo"=mf$yO, "ys"=mf$yS)
print(p[sample(nrow(p), 10),])


### some further tests
# using selection() instead of treatReg()
ms <- selection( yS ~ x + z, yO ~ yS + x, data = dat, 
   mfs = TRUE, mfo = TRUE )
all.equal( m[ ! names( m ) %in% c( "call" , "objectiveFn", "twoStep" ) ], 
   ms[ ! names( ms ) %in% c( "call" , "objectiveFn", "twoStep" ) ] )
# using selection( , type = "treatment" ) instead of treatReg()
mst <- selection( yS ~ x + z, yO ~ yS + x, data = dat, type = "treatment",
   mfs = TRUE, mfo = TRUE )
all.equal( m[ ! names( m ) %in% c( "call" , "objectiveFn", "twoStep" ) ], 
   mst[ ! names( mst ) %in% c( "call" , "objectiveFn", "twoStep" ) ] )
# same treatment variable but a different name in the outcome equation
dat$yD <- dat$yS
try( selection( yS ~ x + z, yO ~ yD + x, data = dat ) )
msD <- selection( yS ~ x + z, yO ~ yD + x, data = dat, type = "treatment",
   mfs = TRUE, mfo = TRUE )
# all.equal( m[ names( m ) != "call" ], msD[ names( msD ) != "call" ] )
# estimate the model as tobit-2 model
try( selection( yS ~ x + z, yO ~ yS + x, data = dat, type = "2" ) )
try( selection( yS ~ x + z, yO ~ yD + x, data = dat, type = "2" ) )

