mvProbitLogLikInternal <- function( yMat, xMat, coef, sigma,
   algorithm, oneSidedGrad, eps, randomSeed, ... ) {

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

   # checking argument 'oneSidedGrad'
   if( length( oneSidedGrad ) != 1 ) {
      stop( "argument 'oneSidedGrad' must be a single logical value" )
   } else if( !is.logical( oneSidedGrad ) ) {
      stop( "argument 'oneSidedGrad' must be logical" )
   }

   # checking argument 'eps'
   if( oneSidedGrad ) {
      if( length( eps ) != 1 ) {
         stop( "argument 'eps' must be a single numeric value" )
      } else if( !is.numeric( eps ) ) {
         stop( "argument 'eps' must be numeric" )
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

   # calculate log likelihood values (for each observation)
   result <- rep( NA, nObs )
   for( i in 1:nObs ){
      ySign <- 2 * yMat[ i, ] - 1
      xBetaTmp <- xBeta[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
      result[ i ] <- log( pmvnormWrap( upper = xBetaTmp, sigma = sigmaTmp, 
         algorithm = algorithm, random.seed = randomSeed, ... ) )
   }

   if( oneSidedGrad ) {
      nDep <- ncol( sigma ) 
      grad <- matrix( NA, nrow = length( result ), 
         ncol = length( coef ) + nDep * ( nDep - 1 ) / 2 )
      for( i in 1:length( coef ) ) {
         coefTmp <- coef
         coefTmp[ i ] <- coef[ i ] + eps
         llTmp <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
            coef = coefTmp, sigma = sigma, algorithm = algorithm,
            oneSidedGrad = FALSE, eps = eps, randomSeed = randomSeed, ... )
         grad[ , i ] <- ( llTmp - result ) / eps
      }
      gradRow <- length( coef )
      for( i in 1:(nDep-1) ) {
         for( j in (i+1):nDep ) {
            gradRow <- gradRow + 1
            sigmaTmp <- sigma
            sigmaTmp[ i, j ] <- sigmaTmp[ j, i ] <- sigma[ j, i ] + eps
            llTmp <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
               coef = coef, sigma = sigmaTmp, algorithm = algorithm,
               oneSidedGrad = FALSE, eps = eps, randomSeed = randomSeed, ... )
            grad[ , gradRow ] <- ( llTmp - result ) / eps
         }
      }
      attr( result, "gradient" ) <- grad
   }

   return( result )
}
