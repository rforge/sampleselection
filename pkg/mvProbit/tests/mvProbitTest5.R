library( "mvProbit" )

## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = rnorm( nObs ),
   x3 = rnorm( nObs ) )

# model coefficients
beta <- cbind( c( -0.4,  0.8, -1.0, -1.4 ),
               c( -0.6, -0.5,  0.6, -1.2 ),
               c(  0.3,  0.9, -1.3,  0.8 ),
               c( -0.3, -1.0,  0.5, -0.8 ),
               c(  1.2, -0.6, -0.7,  1.1 ) )

# covariance matrix of error terms
sigma <- diag( 5 )
sigma[ upper.tri( sigma ) ] <- 
   c( 0.5, 0.4, -0.7, 0.8, -0.5, 0.8, 0, 0.6, -0.2, 0.6 )
sigma <- cov2cor( t( sigma ) %*% sigma )

# all parameters in a vector
allCoef <- c( c( beta ), sigma[ lower.tri( sigma ) ] )

# generate dependent variables
yMatLin <- xMat %*% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:5, sep = "" )
# (yMatLin > 0 )== yMat

# unconditional expectations of dependent variables
yExp <- mvProbitExp( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ) )
print( yExp )
yExpA <- mvProbitExp( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ) )
all.equal( yExp, yExpA )
yExp2 <- pnorm( yMatLin )
all.equal( yExp, as.data.frame( yExp2 ) )


# conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
yExpCond <- mvProbitExp( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE )
print( yExpCond )
yExpCondA <- mvProbitExp( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ), cond = TRUE )
all.equal( yExpCond, yExpCondA )
yExpCond2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ) {
   for( k in 1:ncol( yMat ) ) {
      set.seed( 123 )
      numerator <- pmvnorm( upper = yMatLin[ i, ], sigma = sigma )
      set.seed( 123 )
      denominator <- pmvnorm( upper = yMatLin[ i, -k ], sigma = sigma[ -k, -k ] )
      yExpCond2[ i, k ] <- numerator / denominator
   }
}
all.equal( yExpCond, as.data.frame( yExpCond2 ) )


# conditional expectations of dependent variables
# (assuming that all other dependent variables are as observed)
yExpCondObs <- mvProbitExp( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE )
print( yExpCondObs )
yExpCondObsA <- mvProbitExp( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE )
all.equal( yExpCond, yExpCondA )
yExpCondObs2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ){
   for( k in 1:ncol( yMat ) ) {
      ySign <- 2 * yMat[ i, ] - 1
      ySign[ k ] <- 1
      yLinTmp <- yMatLin[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
      set.seed( 123 )
      numerator <- pmvnorm( upper = yLinTmp, sigma = sigmaTmp )
      set.seed( 123 )
      denominator <- pmvnorm( upper = yLinTmp[ -k ], sigma = sigmaTmp[ -k, -k ] )
      yExpCondObs2[ i, k ] <- numerator / denominator
   }
}
colnames( yExpCondObs2 ) <- names( yExpCondObs )
all.equal( yExpCondObs, as.data.frame( yExpCondObs2 ) )


# calculating log likelihood value(s)
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ) )
print( logLikVal )
logLikValA <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ) )
all.equal( logLikVal, logLikValA )

# calculating log likelihood value(s) with one-sided gradients
logLikValGrad <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   oneSidedGrad = TRUE )
print( logLikValGrad )
logLikValGradA <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ),
   oneSidedGrad = TRUE )
all.equal( logLikValGrad, logLikValGradA )

# calculating log likelihood value(s) with two-sided gradients
llTmp <- function( coef ) {
   betaTmp <- coef[ 1:20 ]
   sigmaTmp <- diag( 5 )
   sigmaTmp[ lower.tri( sigmaTmp ) ] <- coef[ -(1:20) ]
   sigmaTmp[ upper.tri( sigmaTmp ) ] <- t( sigmaTmp )[ upper.tri( sigmaTmp ) ]
   result <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
      coef = betaTmp, sigma = sigmaTmp, 
      data = as.data.frame( cbind( xMat, yMat ) ) )
   return( result )
}
logLikValGrad2 <- numericGradient( llTmp, allCoef )
print( logLikValGrad2 )
attr( logLikValGrad, "gradient" ) / logLikValGrad2 - 1
range( attr( logLikValGrad, "gradient" ) / logLikValGrad2 - 1, na.rm = TRUE )
attr( logLikValGrad, "gradient" ) - logLikValGrad2
range( attr( logLikValGrad, "gradient" ) - logLikValGrad2 )


# calculating marginal effects, unconditional
margEffUnc <- mvProbitMargEff( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ) )
print( margEffUnc )
margEffUncA <- mvProbitMargEff( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ) )
all.equal( margEffUnc, margEffUncA )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are one)
margEffCond <- mvProbitMargEff( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE )
print( margEffCond )
margEffCondA <- mvProbitMargEff( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ), cond = TRUE )
all.equal( margEffCond, margEffCondA )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are as observed)
margEffCondObs <- mvProbitMargEff( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE )
print( margEffCondObs )
margEffCondObsA <- mvProbitMargEff( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ), cond = TRUE )
all.equal( margEffCondObs, margEffCondObsA )

