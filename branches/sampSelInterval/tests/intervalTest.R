library( "sampleSelection" )
library( "mvtnorm" )
options( digits = 2 )

nObs <- 300
betaS <- c( 1, 1, -1 )
betaO <- c( 10, 4 )
rho <- 0.4
sigma <- 5
bound <- c(-Inf,5,15,Inf)


set.seed(123)
dat <- data.frame( x1 = rnorm( nObs ), x2 = rnorm( nObs ) )
vcovMat <- matrix( c( 1, rho*sigma, rho*sigma, sigma^2 ), nrow = 2 )
eps <- rmvnorm( nObs, sigma = vcovMat )
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

start <- c( betaS, betaO, log( sigma ), atan( rho ) )
   # the correct starting value of logSigmaSq would be: log( sigma^2 )
names( start ) <- c( "betaS0", "betaS1", "betaS2", "betaO0", "betaO2",
   "logSigmaSq", "atanRho" )

res <- sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = bound, 
    start = start, printLevel = 1 )

print( res )
print( round( coef( res ), 2 ) )
print( round( coef( summary( res ) ), 2 ) )
print( res$start )

# add derived coefficients
coefAll <- c( coef( res ),
   sigma = unname( sqrt( exp( coef( res )[ "logSigmaSq" ] ) ) ),
   sigmaSq = unname( exp( coef( res )[ "logSigmaSq" ] ) ),
   rho = unname( tan( coef( res )[ "atanRho" ] ) ) )
print( round( coefAll, 2 ) )

# jacobian
jac <- cbind( diag( length( coef( res ) ) ),
   matrix( 0, length( coef( res ) ), 3 ) )
rownames( jac ) <- names( coef( res ) )
colnames( jac ) <- c( names( coef( res ) ), "sigma", "sigmaSq", "rho" )
jac[ "logSigmaSq", "sigma" ] <- sqrt( exp( coef( res )[ "logSigmaSq" ] ) ) / 2
jac[ "logSigmaSq", "sigmaSq" ] <- exp( coef( res )[ "logSigmaSq" ] )
jac[ "atanRho", "rho" ] <- 1 + ( tan( coef( res )[ "atanRho" ] ) )^2
vcovAll <- t( jac ) %*% vcov( res ) %*% jac
print( round( vcovAll, 2 ) )
print( round( cov2cor( vcovAll ), 2 ) )
print( coefTable( coefAll, sqrt( diag( vcovAll ) ) ) ) 

maxLik:::summary.maxLik( res )

# function that returns log-likelihood values
intLogLik <- function( param ) {
   ll <- sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = bound, 
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


# tests with automatically generated starting values (ML estimation)
resMl <- sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = bound, 
   start = "ml", printLevel = 1 )
print( resMl )
print( round( coef( resMl ), 2 ) )
print( round( coef( summary( resMl ) ), 2 ) )
print( resMl$start )


# tests with automatically generated starting values (2-step estimation)
res2s <- sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = bound, 
   start = "2step", printLevel = 1 )
print( res2s )
print( round( coef( res2s ), 2 ) )
print( round( coef( summary( res2s ) ), 2 ) )
print( res2s$start )


# tests with incorrectly specified starting values
try( sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = bound, 
   start = "wrong" ) )
try( sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = bound, 
   start = rep( 1, 11 ) ) )

# tests with incorrectly specified boundaries
try( sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = 1:6, 
   start = start ) )
try( sampleSelection:::tobit2Intfit( YS, XS, YO, XO, boundaries = 4:1, 
   start = start ) )
