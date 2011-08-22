mvProbitExpInternal <- function( yMat, xMat, coef, sigma,
   cond = FALSE, random.seed = 123, ... ) {

   # checking argument 'cond'
   if( !is.logical( cond ) ) {
      stop( "argument 'cond' must be logical" )
   } else if( length( cond ) != 1 ) {
      stop( "argument 'cond' must be a single logical values" )
   }

   # checking argument 'random.seed'
   if( !is.numeric( random.seed ) ) {
      stop( "argument 'random.seed' must be numerical" )
   } else if( length( random.seed ) != 1 ) {
      stop( "argument 'random.seed' must be a single numerical values" )
   }

   # checking argument 'sigma'
   if( !is.matrix( sigma ) ) {
      stop( "argument 'sigma' must be a matrix" )
   } else if( nrow( sigma ) != ncol( sigma ) ) {
      stop( "argument 'sigma' must be a quadratic matrix" )
   } else if( !isSymmetric( sigma ) ) {
      stop( "argument 'sigma' must be a symmetric matrix" )
   } else if( !is.null( yMat ) ) {
      if( ncol( sigma ) != ncol( yMat ) ) {
         stop( "the number of dependent variables specified in argument",
            " 'formula' must be equal to the number of rows and colums",
            " of the matrix specified by argument 'sigma'" )
      }
   }

   # number of dependent variables
   nDep <- ncol( sigma )

   # number of regressors
   nReg <- ncol( xMat )

   # number of coefficients
   nCoef <- nDep * nReg

   # number of observations
   nObs <- nrow( xMat )

   # checking argument 'coef'
   if( !is.vector( coef, mode = "numeric" ) ) {
      stop( "argument 'coef' must be a numeric vector" )
   } else if( length( coef ) != nCoef ) {
      stop( "argument coef must have ", nCoef, " elements" )
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

   # save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by pmvnorm)
   set.seed( random.seed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   if( cond ) {
      # conditional expectations
      result <- matrix( NA, nrow = nObs, ncol = nDep )
      if( is.null( yMat ) ) {
         # assuming that all other dependent variables are one
         for( i in 1:nObs ) {
            for( k in 1:nDep ) {
               result[ i, k ] <- 
                  pmvnorm( upper = xBeta[ i, ], sigma = sigma, ... ) / 
                  pmvnorm( upper = xBeta[ i, -k ], sigma = sigma[ -k, -k ], ... )
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
                  pmvnorm( upper = xBetaTmp, sigma = sigmaTmp, ... ) / 
                  pmvnorm( upper = xBetaTmp[ -k ], sigma = sigmaTmp[ -k, -k ], ... )
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
