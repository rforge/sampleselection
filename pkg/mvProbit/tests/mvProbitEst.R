library( "mvProbit" )
library( "miscTools" )

## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 50

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = rnorm( nObs ) )

# coefficients
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
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5 )
summary( estResultBHHH )

# estimation with the BHHH algorithm, one-sided gradients
estResultBHHH1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5,
   oneSidedGrad = TRUE )
summary( estResultBHHH1 )
all.equal( estResultBHHH, estResultBHHH1 )

# estimation with the BFGS algorithm, two-sided gradients
estResultBFGS <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5 )
summary( estResultBFGS )

# estimation with the BFGS algorithm, one-sided gradients
estResultBFGS1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, oneSidedGrad = TRUE )
summary( estResultBFGS1 )
all.equal( estResultBFGS, estResultBFGS1 )

# estimation with the BFGS algorithm, one-sided gradients, no starting values
estResultBFGS1a <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, oneSidedGrad = TRUE )
summary( estResultBFGS1a )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for beta
estResultBFGS1b <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, oneSidedGrad = TRUE )
summary( estResultBFGS1b )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for sigma
estResultBFGS1s <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, oneSidedGrad = TRUE )
summary( estResultBFGS1s )

# estimation with the BFGS algorithm, Miwa algorithm for obtaining integrals
estResultBFGSm <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, algorithm = Miwa( steps = 64 ) )
summary( estResultBFGSm )
all.equal( estResultBFGS, estResultBFGSm )

# estimation with the BFGS algorithm, GHK algorithm for obtaining integrals
estResultBFGSg <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, algorithm = "GHK" )
summary( estResultBFGSg )
all.equal( estResultBFGS, estResultBFGSg )
all.equal( estResultBFGSm, estResultBFGSg )

# estimation with the Nelder-Mead algorithm
estResultNM <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "NM", reltol = 0.05 )
summary( estResultNM )

