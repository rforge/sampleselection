library( "mvProbit" )
library( "miscTools" )

## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = as.numeric( rnorm( nObs ) > 0 ),
   x3 = rnorm( nObs ),
   x4 = rnorm( nObs ) )

# coefficients
beta <- cbind( c(  0.8,  1.2, -1.0,  1.4, -0.8 ),
               c( -0.6,  1.0,  0.6, -1.2, -1.6 ),
               c(  0.5, -0.6, -0.7,  1.1,  1.2 ) )

# covariance matrix of error terms
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat %*% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )
# (yMatLin > 0 )== yMat

# unconditional expectations of dependent variables
yExp <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ) )
print( yExp )
yExp2 <- pnorm( yMatLin )
all.equal( yExp, as.data.frame( yExp2 ) )

# conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
yExpCond <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE )
print( yExpCond )
set.seed( 123 )
yExpCond2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ) {
   for( k in 1:ncol( yMat ) ) {
      yExpCond2[ i, k ] <- pmvnorm( upper = yMatLin[ i, ], sigma = sigma ) / 
         pmvnorm( upper = yMatLin[ i, -k ], sigma = sigma[ -k, -k ] )
   }
}
all.equal( yExpCond, as.data.frame( yExpCond2 ) )

# conditional expectations of dependent variables
# (assuming that all other dependent variables are as observed)
yExpCondObs <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE )
print( yExpCondObs )
set.seed( 123 )
yExpCondObs2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ){
   for( k in 1:ncol( yMat ) ) {
      ySign <- 2 * yMat[ i, ] - 1
      ySign[ k ] <- 1
      yLinTmp <- yMatLin[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
      yExpCondObs2[ i, k ] <- pmvnorm( upper = yLinTmp, sigma = sigmaTmp ) / 
         pmvnorm( upper = yLinTmp[ -k ], sigma = sigmaTmp[ -k, -k ] )
   }
}
all.equal( yExpCondObs, as.data.frame( yExpCondObs2 ) )

# unconditional expectations of dependent variables by simulation
nSim <- 10000
ySim <- array( NA, c( nObs, ncol( yMat ), nSim ) )
for( s in 1:nSim ) {
   ySim[ , , s ] <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
}
yExpSim <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ) {
   yExpSim[ i, ] <- rowSums( ySim[ i, , ] ) / nSim
}
print( yExpSim )
print( yExpSim - as.matrix( yExp ) )

# calculating log likelihood value(s)
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ) )
print( logLikVal )

# calculating marginal effects, unconditional
margEffUnc <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ) )
print( margEffUnc )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are one)
margEffCond <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE )
print( margEffCond )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are as observed)
margEffCondObs <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE )
print( margEffCondObs )

