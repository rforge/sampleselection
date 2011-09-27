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
# now with explicitly specifying the algorithm
yExpCond3 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz )
all.equal( yExpCond, yExpCond3 )
identical( yExpCond, yExpCond3 )
# now with integrals obtained by the Miwa algorithm
yExpCond4 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = Miwa )
all.equal( yExpCond, yExpCond4 )
# now with integrals obtained by the Miwa algorithm, less precise
yExpCond5 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = Miwa( steps = 32 ) )
all.equal( yExpCond4, yExpCond5 )
# now with integrals obtained by the TVPACK algorithm
yExpCond6 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = TVPACK )
all.equal( yExpCond, yExpCond6 )
# now with integrals obtained by the TVPACK algorithm, less precise
yExpCond7 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = TVPACK( abseps = 0.5 ) )
all.equal( yExpCond6, yExpCond7 )


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
# now with explicitly specifying the algorithm
yExpCondObs3 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = GenzBretz )
all.equal( yExpCondObs, yExpCondObs3 )
identical( yExpCondObs, yExpCondObs3 )
# now with integrals obtained by the Miwa algorithm
yExpCondObs4 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = Miwa )
all.equal( yExpCondObs, yExpCondObs4 )
# now with integrals obtained by the Miwa algorithm, less precise
yExpCondObs5 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = Miwa( steps = 32 ) )
all.equal( yExpCondObs4, yExpCondObs5 )
# now with integrals obtained by the TVPACK algorithm
yExpCondObs6 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = TVPACK )
all.equal( yExpCondObs, yExpCondObs6 )
# now with integrals obtained by the TVPACK algorithm, less precise
yExpCondObs7 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = TVPACK( abseps = 0.5 ) )
all.equal( yExpCondObs6, yExpCondObs7 )


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

# for testing state of random number generator
rnorm( 4 )

# calculating log likelihood value(s)
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ) )
print( logLikVal )

# calculating log likelihood value(s) with one-sided gradients
logLikValGrad <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   oneSidedGrad = TRUE )
print( logLikValGrad )

# calculating log likelihood value(s) with two-sided gradients
llTmp <- function( coef ) {
   betaTmp <- coef[ 1:15 ]
   sigmaTmp <- diag( 3 )
   sigmaTmp[ upper.tri( sigmaTmp ) ] <- coef[ -(1:15) ]
   sigmaTmp[ lower.tri( sigmaTmp ) ] <- t( sigmaTmp )[ lower.tri( sigmaTmp ) ]
   result <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
      coef = betaTmp, sigma = sigmaTmp, 
      data = as.data.frame( cbind( xMat, yMat ) ) )
   return( result )
}
allCoef <- c( c( beta ), sigma[ upper.tri( sigma ) ] )
logLikValGrad2 <- numericGradient( llTmp, allCoef )
print( logLikValGrad2 )
attr( logLikValGrad, "gradient" ) / logLikValGrad2 - 1
range( attr( logLikValGrad, "gradient" ) / logLikValGrad2 - 1, na.rm = TRUE )
attr( logLikValGrad, "gradient" ) - logLikValGrad2
range( attr( logLikValGrad, "gradient" ) - logLikValGrad2 )

# for testing state of random number generator
rnorm( 4 )

# calculating marginal effects, unconditional
margEffUnc <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ) )
print( margEffUnc )

# for testing state of random number generator
rnorm( 4 )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are one)
margEffCond <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE )
print( margEffCond )
# now with integrals obtained by the Miwa algorithm, reduced precision
margEffCond1 <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = Miwa( steps = 32 ) )
print( margEffCond1 )
all.equal( margEffCond, margEffCond1 )

# for testing state of random number generator
rnorm( 4 )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are as observed)
margEffCondObs <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE )
print( margEffCondObs )
# now with integrals obtained by the Miwa algorithm, reduced precision
margEffCondObs1 <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = Miwa( steps = 32 ) )
print( margEffCondObs1 )
all.equal( margEffCondObs, margEffCondObs1 )

# for testing state of random number generator
rnorm( 4 )

