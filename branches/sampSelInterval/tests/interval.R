library( "sampleSelection" )
library( "mvtnorm" )
options( digits = 2 )

nObs <- 300
betaS <- c( 1, 1, -1 )
betaO <- c( 10, 4 )
rho <- 0.4
sigma2 <- 5
bound <- c(-Inf,5,15,Inf)


set.seed(123)
dat <- data.frame( x1 = rnorm( nObs ), x2 = rnorm( nObs ) )
Sigma <- matrix( c( 1, rho*sigma2, rho*sigma2, sigma2^2 ), nrow = 2 )
eps <- rmvnorm( nObs, sigma = Sigma )
dat$epsS <- eps[,1]
dat$epsO <- eps[,2]
dat$yS <- with( dat, betaS[1] + betaS[2] * x1 + betaS[3] * x2 + epsS ) > 0
# table(dat$yS)
dat$yOu <- with( dat, betaO[1] + betaO[2] * x1 + epsO )
dat$yOu[ !dat$yS ] <- NA
hist( dat$yOu )
dat$yO <- cut( dat$yOu, bound )

YS <- dat$yS
XS <- cbind( 1, dat$x1, dat$x2 )
YO <- as.numeric( dat$yO )
XO <- cbind( 1, dat$x1 )

start <- c( betaS, betaO, rho, sigma2 )
names( start ) <- c( "betaS0", "betaS1", "betaS2", "betaO0", "betaO2",
   "atanrho", "logsigmaSq" )

res <- sampleSelection:::intervalfit( YS, XS, YO, XO, boundaries = bound, 
    start = start, printLevel = 1 )

print( res )
print( round( coef( res ), 2 ) )
print( round( coef( summary( res ) ), 2 ) )
maxLik:::summary.maxLik( res )

# function that returns log-likelihood values
intLogLik <- function( param ) {
   ll <- sampleSelection:::intervalfit( YS, XS, YO, XO, boundaries = bound, 
      start = param, returnLogLikStart = TRUE )
   return( ll )
}


# log-likelihood values and their gradients at starting values of parameters
logLikStart <- intLogLik( start )
print( c( logLikStart ) )
gradStart <- attr( logLikStart, "gradient" )
colnames( gradStart ) <- names( start )
print( gradStart )
# numeric gradients
gradStartNum <- numericGradient( intLogLik, t0 = start )
colnames( gradStartNum ) <- names( start )
all.equal( as.data.frame( gradStart ), as.data.frame( gradStartNum ) )
library( "miscTools" )
for( i in 1:ncol( gradStart ) ) {
   compPlot( gradStart[ , i ], gradStartNum[ , i ], main = names( start )[i],
      col = ifelse( !YS, "black",
         ifelse( YO == 1, "blue",
            ifelse( YO == 2, "blueviolet",
               ifelse( YO == 3, "red", "green" ) ) ) ) )
}

# log-likelihood values and their gradients at estimated parameters
logLikEst <- intLogLik( coef( res ) )
print( c( logLikEst ) )
gradEst <- gradStart <- attr( logLikEst, "gradient" )
colnames( gradEst ) <- names( start )
print( gradEst )


