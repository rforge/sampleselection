mvProbitLogLikInternal <- function( yMat, xMat, coef, sigma,
   algorithm, nGHK, oneSidedGrad, eps, randomSeed, ... ) {

   # number of regressors
   nReg <- ncol( xMat )

   # checking and preparing model coefficients and correlation coefficients
   coef <- mvProbitPrepareCoef( yMat = yMat, nReg = nReg, coef = coef, 
      sigma = sigma )

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
   nDep <- ncol( coef$sigma )

   # number of observations
   nObs <- nrow( xMat )

   # calculating linear predictors
   xBeta <- matrix( NA, nrow = nObs, ncol = nDep )
   for( i in 1:nDep ) {
      xBeta[ , i ] <- xMat %*% coef$betaEq[[ i ]]
   }

   # calculate log likelihood values (for each observation)
   result <- rep( NA, nObs )
   for( i in 1:nObs ){
      ySign <- 2 * yMat[ i, ] - 1
      xBetaTmp <- xBeta[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% coef$sigma %*% diag( ySign )
      result[ i ] <- log( pmvnormWrap( upper = xBetaTmp, sigma = sigmaTmp, 
         algorithm = algorithm, nGHK = nGHK, random.seed = randomSeed, ... ) )
   }

   if( oneSidedGrad ) {
      allCoef <- c( coef$beta, coef$sigma[ lower.tri( coef$sigma ) ] )
      grad <- matrix( NA, nrow = length( result ), 
         ncol = length( allCoef ) )
      for( i in 1:nDep ) {
         # gradients of intercepts
         coefTmp <- allCoef
         InterceptNo <- ( i - 1 ) * nReg + 1
         coefTmp[ InterceptNo ] <- allCoef[ InterceptNo ] + eps
         llTmp <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
            coef = coefTmp, sigma = NULL, algorithm = algorithm, nGHK = nGHK,
            oneSidedGrad = FALSE, eps = eps, randomSeed = randomSeed, ... )
         grad[ , InterceptNo ] <- ( llTmp - result ) / eps
         # gradients of other variables
         if( nReg > 1 ) {
            for( j in 2:nReg ) {
               grad[ , InterceptNo + j - 1 ] <- 
                  grad[ , InterceptNo ] * xMat[ , j ] 
            }
         }
      }
      for( i in ( nDep * nReg + 1 ):length( allCoef ) ) {
         coefTmp <- allCoef
         coefTmp[ i ] <- allCoef[ i ] + eps
         llTmp <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
            coef = coefTmp, sigma = NULL, algorithm = algorithm, nGHK = nGHK,
            oneSidedGrad = FALSE, eps = eps, randomSeed = randomSeed, ... )
         grad[ , i ] <- ( llTmp - result ) / eps
      }
      colnames( grad ) <- mvProbitCoefNames( nDep = nDep, nReg = nReg )
      attr( result, "gradient" ) <- grad
   }

   return( result )
}
