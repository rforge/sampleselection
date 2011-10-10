mvProbitLogLikInternal <- function( yMat, xMat, coef, sigma,
   algorithm, nGHK, oneSidedGrad, eps, randomSeed, ... ) {

   # number of regressors
   nReg <- ncol( xMat )

   # checking and preparing coefficients and correlation coefficients
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
      nDep <- ncol( coef$sigma ) 
      grad <- matrix( NA, nrow = length( result ), 
         ncol = length( coef$beta ) + nDep * ( nDep - 1 ) / 2 )
      for( i in 1:length( coef$beta ) ) {
         coefTmp <- coef$beta
         coefTmp[ i ] <- coef$beta[ i ] + eps
         llTmp <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
            coef = coefTmp, sigma = coef$sigma, algorithm = algorithm, nGHK = nGHK,
            oneSidedGrad = FALSE, eps = eps, randomSeed = randomSeed, ... )
         grad[ , i ] <- ( llTmp - result ) / eps
      }
      gradRow <- length( coef$beta )
      for( i in 1:(nDep-1) ) {
         for( j in (i+1):nDep ) {
            gradRow <- gradRow + 1
            sigmaTmp <- coef$sigma
            sigmaTmp[ i, j ] <- sigmaTmp[ j, i ] <- coef$sigma[ j, i ] + eps
            llTmp <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
               coef = coef$beta, sigma = sigmaTmp, algorithm = algorithm, nGHK = nGHK,
               oneSidedGrad = FALSE, eps = eps, randomSeed = randomSeed, ... )
            grad[ , gradRow ] <- ( llTmp - result ) / eps
         }
      }
      colnames( grad ) <- mvProbitCoefNames( nDep = nDep, nReg = nReg )
      attr( result, "gradient" ) <- grad
   }

   return( result )
}
