mvProbitLogLikInternal <- function( yMat, xMat, coef, sigma,
   randomSeed = 123, ... ) {

   # checking argument 'random.seed' / 'randomSeed'
   if( !is.numeric( randomSeed ) ) {
      stop( "argument 'random.seed' must be numerical" )
   } else if( length( randomSeed ) != 1 ) {
      stop( "argument 'random.seed' must be a single numerical values" )
   }

   # checking argument 'sigma'
   if( !is.matrix( sigma ) ) {
      stop( "argument 'sigma' must be a matrix" )
   } else if( nrow( sigma ) != ncol( sigma ) ) {
      stop( "argument 'sigma' must be a quadratic matrix" )
   } else if( !isSymmetric( sigma ) ) {
      stop( "argument 'sigma' must be a symmetric matrix" )
   } else if( any( abs( diag( sigma ) - 1 ) > 1e-7 ) ) {
      stop( "argument 'sigma' must have ones on its diagonal" )
   } else if( ncol( sigma ) != ncol( yMat ) ) {
      stop( "the number of dependent variables specified in argument",
         " 'formula' must be equal to the number of rows and colums",
         " of the matrix specified by argument 'sigma'" )
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
   set.seed( randomSeed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   # calculate log likelihood values (for each observation)
   result <- rep( NA, nObs )
   for( i in 1:nObs ){
      ySign <- 2 * yMat[ i, ] - 1
      xBetaTmp <- xBeta[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
      result[ i ] <- log( pmvnorm( upper = xBetaTmp, sigma = sigmaTmp, ... ) )
   }

   return( result )
}
