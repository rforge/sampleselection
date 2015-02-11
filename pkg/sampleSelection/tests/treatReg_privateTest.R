### Testing treatreg.
### These are testst that are not supposed to be included in CRAN

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
set.seed(1)
options(digits=4)
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