mvProbit <- function( formula, data, start = NULL, startSigma = NULL,
   method = "BFGS", finalHessian = "BHHH", algorithm = GenzBretz(), nGHK = 1000,
   oneSidedGrad = FALSE, eps = 1e-6, random.seed = 123, ... ) {

   # checking argument 'formula'
   if( is.list( formula ) ) {
      stop( "using different regressors for the dependent variables",
         " has not been implemented yet. Sorry!" )
   } else if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   }

   # checking argument 'data'
   if( !is.data.frame( data ) ) {
      stop( "argument 'data' must be a data frame" )
   }

   # preparing model matrix
   mc <- match.call( expand.dots = FALSE )
   m <- match( "data", names( mc ), 0 )
   mf <- mc[ c( 1, m ) ]
   mf$formula <- formula
   attributes( mf$formula ) <- NULL
   mf$na.action <- na.pass
   mf[[ 1 ]] <- as.name( "model.frame" )
   mf <- eval( mf, parent.frame() )
   mt <- attr( mf, "terms" )
   xMat <- model.matrix( mt, mf )

   # preparing model response
   yMat <- model.response( mf )
   if( !is.matrix( yMat ) ) {
      stop( "at least two dependent variables",
         " must be specified in argument 'formula'",
         " (e.g. by 'cbind( y1, y2 ) ~ ...')" )
   } else if( !all( yMat %in% c( 0, 1,TRUE, FALSE ) ) ) {
      stop( "all dependent variables must be either 0, 1, TRUE, or FALSE" )
   }

   # number of dependent variables
   nDep <- ncol( yMat )

   # number of regressors
   nReg <- ncol( xMat )

   # number of model coefficients
   nCoef <- nDep * nReg

   # number of observations
   nObs <- nrow( xMat )

   # obtaining starting values for model coefficients if they are not specified
   if( is.null( start ) ) {
      uvProbit <- list()
      for( i in 1:nDep ) {
         uvProbit[[ i ]] <- glm( yMat[ , i ] ~ xMat - 1, 
            family = binomial( link = "probit" ) )
         start <- c( start, coef( uvProbit[[ i ]] ) )
      }
   }

   # obtaining starting values for correlations if they are not specified
   if( is.null( startSigma ) && length( start ) != nCoef + nDep * ( nDep - 1 ) / 2 ) {
      yHat <- matrix( NA, nrow = nObs, ncol = nDep )
      for( i in 1:nDep ) {
         yHat[ , i ] <- pnorm( xMat %*% 
            start[ ( ( i - 1 ) * nReg + 1 ):( i * nReg ) ] )
      }
      startSigma <- cor( yMat - yHat )
   }

   # checking and preparing model coefficients and correlation coefficients
   coef <- mvProbitPrepareCoef( yMat = yMat, nReg = nReg, coef = start, 
      sigma = startSigma )

   # starting values
   start <- c( coef$beta, coef$sigma[ lower.tri( coef$sigma ) ] )
   names( start ) <- mvProbitCoefNames( nDep = nDep, nReg = nReg )

   # wrapper function for maxLik for calling mvProbitLogLikInternal
   logLik <- function( param, yMat, xMat,
      llAlgorithm, llNGHK, llOneSidedGrad, llEps, llRandom.seed, ... ) {

      logLikVal <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
         coef = param, sigma = NULL, algorithm = llAlgorithm, nGHK = llNGHK,
         oneSidedGrad = llOneSidedGrad, eps = llEps, randomSeed = llRandom.seed, 
         ... )

      return( logLikVal )
   }

   result <- maxLik( logLik = logLik, start = start, method = method, 
      finalHessian = finalHessian, 
      yMat = yMat, xMat = xMat, llAlgorithm = algorithm, llNGHK = nGHK,
      llOneSidedGrad = oneSidedGrad, llEps = eps, llRandom.seed = random.seed,
      ... )

   # return also some other useful information
   result$call <- match.call()
   result$start <- start
   result$nDep <- nDep
   result$nReg <- nReg
   result$nObs <- nObs

   class( result ) <- c( "mvProbit", class( result ) )

   return( result )
}
