library( "mvProbit" )

## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 50

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = rnorm( nObs ) )

# model coefficients
beta <- cbind( c(  0.8,  1.2, -0.8 ),
               c( -0.6,  1.0, -1.6 ),
               c(  0.5, -0.6,  1.2 ) )

# covariance matrix of error terms
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat %*% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )

# estimation with the BHHH algorithm, two-sided gradients
estResultBHHH <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, method = "BHHH",
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5 )
print( estResultBHHH )
summary( estResultBHHH )
estResultBHHHA <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta, sigma[ lower.tri( sigma ) ] ), method = "BHHH",
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5 )
all.equal( estResultBHHH, estResultBHHHA )

# estimation with the BHHH algorithm, one-sided gradients
estResultBHHH1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, method = "BHHH",
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5,
   oneSidedGrad = TRUE )
print( estResultBHHH1 )
summary( estResultBHHH1 )
all.equal( estResultBHHH, estResultBHHH1 )

# estimation with the BFGS algorithm, two-sided gradients
estResultBFGS <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5 )
print( estResultBFGS )
summary( estResultBFGS )

# estimation with the BFGS algorithm, one-sided gradients
estResultBFGS1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, oneSidedGrad = TRUE )
print( estResultBFGS1 )
summary( estResultBFGS1 )
all.equal( estResultBFGS, estResultBFGS1 )

# estimation with the BFGS algorithm, one-sided gradients, no starting values
estResultBFGS1a <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, oneSidedGrad = TRUE )
print( estResultBFGS1a )
summary( estResultBFGS1a )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for beta
estResultBFGS1b <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   startSigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, oneSidedGrad = TRUE )
print( estResultBFGS1b )
summary( estResultBFGS1b )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for sigma
estResultBFGS1s <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, oneSidedGrad = TRUE )
print( estResultBFGS1s )
summary( estResultBFGS1s )

# estimation with the BFGS algorithm, Miwa algorithm for obtaining integrals
estResultBFGSm <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, algorithm = Miwa( steps = 64 ) )
print( estResultBFGSm )
summary( estResultBFGSm )
all.equal( estResultBFGS, estResultBFGSm )

# estimation with the BFGS algorithm, GHK algorithm for obtaining integrals
estResultBFGSg <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, algorithm = "GHK" )
print( estResultBFGSg )
summary( estResultBFGSg )
all.equal( estResultBFGS, estResultBFGSg )
all.equal( estResultBFGSm, estResultBFGSg )

# estimation with the Nelder-Mead algorithm
estResultNM <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "NM", reltol = 0.05 )
print( estResultNM )
summary( estResultNM )

# marginal effects based on estimated coefficients with covariance matrix
# unconditional marginal effects
margEffUnc <- margEff( estResultBFGS, calcVCov = TRUE )
print( margEffUnc )
print( attr( margEffUnc, "vcov" )[ 1:5, , ] )
print( drop( attr( margEffUnc, "vcov" )[ nObs, , ] ) )
summary( margEffUnc )

# conditional marginal effects
# (assuming that all other dependent variables are as observed)
margEffCondObs <- margEff( estResultBFGS, cond = TRUE )
print( margEffCondObs )

# conditional marginal effects with covariance matrix at sample mean
# (assuming that all other dependent variables are at there modal values)
margEffCondObsCov <- margEff( estResultBFGS, cond = TRUE,
   data = as.data.frame( t( c( colMedians( yMat * 1 ), colMeans( xMat ) ) ) ), 
   calcVCov = TRUE )
print( margEffCondObsCov )
print( attr( margEffCondObsCov, "vcov" ) )
print( drop( attr( margEffCondObsCov, "vcov" ) ) )
summary( margEffCondObsCov )

# conditional marginal effects
# (assuming that all other dependent variables are one)
margEffCondOne <- margEff( estResultBFGS, cond = TRUE, othDepOne = TRUE )
print( margEffCondOne )

# conditional marginal effects with covariance matrix at sample mean
# (assuming that all other dependent variables are one)
margEffCondOneCov <- margEff( estResultBFGS, cond = TRUE, othDepOne = TRUE,
   data = as.data.frame( t( colMeans( xMat ) ) ), calcVCov = TRUE )
print( margEffCondOneCov )
print( attr( margEffCondOneCov, "vcov" ) )
print( drop( attr( margEffCondOneCov, "vcov" ) ) )
summary( margEffCondOneCov )


