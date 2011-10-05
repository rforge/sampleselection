mvProbitExpInternal <- function( yMat, xMat, coef, sigma,
   cond, algorithm, nGHK, random.seed, ... ) {

   # checking argument 'cond'
   if( !is.logical( cond ) ) {
      stop( "argument 'cond' must be logical" )
   } else if( length( cond ) != 1 ) {
      stop( "argument 'cond' must be a single logical values" )
   }

   # number of regressors
   nReg <- ncol( xMat )

   # number of observations
   nObs <- nrow( xMat )

   # checking argument 'coef'
   if( !is.vector( coef, mode = "numeric" ) ) {
      stop( "argument 'coef' must be a numeric vector" )
   }

   if( !is.null( sigma ) ) {
      # checking argument 'sigma'
      if( !is.matrix( sigma ) ) {
         stop( "argument 'sigma' must be a matrix" )
      } else if( nrow( sigma ) != ncol( sigma ) ) {
         stop( "argument 'sigma' must be a quadratic matrix" )
      } else if( !isSymmetric( sigma ) ) {
         stop( "argument 'sigma' must be a symmetric matrix" )
      } else if( any( abs( diag( sigma ) - 1 ) > 1e-7 ) ) {
         stop( "argument 'sigma' must have ones on its diagonal" )
      } else if( !is.null( yMat ) ) {
         if( ncol( sigma ) != ncol( yMat ) ) {
            stop( "the number of dependent variables specified in argument",
               " 'formula' must be equal to the number of rows and colums",
               " of the matrix specified by argument 'sigma'" )
         }
      }
      # number of dependent variables
      nDep <- ncol( sigma )
      # number of coefficients
      nCoef <- nDep * nReg
      if( length( coef ) != nCoef ) {
         stop( "given that argument 'sigma' has been specified",
            " argument coef must have ", nCoef, " elements" )
      }
   } else {
      if( !is.null( yMat ) ) {
         # number of dependent variables
         nDep <- ncol( yMat )
         # number of coefficients
         nCoef <- nDep * nReg
         # number of parameters including sigma
         nCoefSigma <- nCoef + nDep * ( nDep - 1 ) / 2
         if( length( coef ) != nCoefSigma ) {
            stop( "given that argument 'sigma' is 'NULL'",
               " argument coef must have ", nCoefSigma, " elements" )
         }
      } else {
         # number of dependent variables
         nDep <- round( - nReg + 0.5 + 
            sqrt( ( nReg - 0.5 )^2 + 2 * length( coef ) ) )
         # number of coefficients
         nCoef <- nDep * nReg
         # number of parameters including sigma
         nCoefSigma <- nCoef + nDep * ( nDep - 1 ) / 2
         if( length( coef ) != nCoefSigma ) {
            stop( "given that argument 'sigma' is 'NULL'",
               " argument coef must have ", nCoefSigma, " elements",
               " if the model has ", nDep, " dependent variables" )
         }
      }
      # extracting correlation coefficients from 'coef' if they are there
      sigma <- diag( nDep )
      sigma[ lower.tri( sigma ) ] <- coef[ -c( 1:nCoef ) ]
      sigma[ upper.tri( sigma ) ] <- t( sigma )[ upper.tri( sigma ) ]
      coef <- coef[ 1:nCoef ]
   }

   # separating coefficients for different equations
   betaEq <- list()
   for( i in 1:nDep ) {
      betaEq[[ i ]] <- coef[ ( ( i - 1 ) * nReg + 1 ):( i * nReg ) ]
   }

   # calculating linear predictors
   xBeta <- matrix( NA, nrow = nObs, ncol = nDep )
   for( i in 1:nDep ) {
      xBeta[ , i ] <- xMat %*% betaEq[[ i ]]
   }

   if( cond ) {
      # conditional expectations
      result <- matrix( NA, nrow = nObs, ncol = nDep )
      if( is.null( yMat ) ) {
         # assuming that all other dependent variables are one
         for( i in 1:nObs ) {
            for( k in 1:nDep ) {
               result[ i, k ] <- 
                  pmvnormWrap( upper = xBeta[ i, ], sigma = sigma,
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... ) /
                  pmvnormWrap( upper = xBeta[ i, -k ], sigma = sigma[ -k, -k ],
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... )
            }
         }
      } else {
         # assuming that all other dependent variables are as observed
         for( i in 1:nObs ){
            for( k in 1:nDep ) {
               ySign <- 2 * yMat[ i, ] - 1
               ySign[ k ] <- 1
               xBetaTmp <- xBeta[ i, ] * ySign
               sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
               result[ i, k ] <- 
                  pmvnormWrap( upper = xBetaTmp, sigma = sigmaTmp,
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... ) / 
                  pmvnormWrap( upper = xBetaTmp[ -k ], sigma = sigmaTmp[ -k, -k ],
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... )
            }
         }
      }
   } else {
      result <- pnorm( xBeta )
   }

   if( !is.null( yMat ) ) {
      colnames( result ) <- colnames( yMat )
   }

   result <- as.data.frame( result )

   return( result )
}
