mvProbit <- function( formula, coef, sigma, data,
   method = "BHHH", finalHessian = "BHHH",
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

   # starting values
   start <- c( coef, sigma[ upper.tri( sigma ) ] )
   names( start ) <- c( 
      paste( "b", rep( 1:nDep, each = nReg ), rep( 0:(nReg-1), nDep ), 
         sep = "_" ),
      matrix( paste( "R", rep( 1:nDep, each = nDep ), rep( 1:nDep, nDep ),
         sep = "_" ), nrow = nDep, byrow = TRUE )[ upper.tri( diag(nDep) ) ] )

   # wrapper function for maxLik for calling mvProbitLogLikInternal
   logLik <- function( param, yMat, xMat, randomSeed, nCoef, nDep, 
      llOneSidedGrad, llEps, ... ) {

      coef <- param[ 1:nCoef ]
      sigma <- diag( nDep )
      sigma[ upper.tri( sigma ) ] <- param[ -(1:nCoef) ]
      sigma[ lower.tri( sigma ) ] <- t( sigma )[ lower.tri( sigma ) ]
      
      logLikVal <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
         coef = coef, sigma = sigma, randomSeed = randomSeed, 
         oneSidedGrad = llOneSidedGrad, eps = llEps, ... )

      return( logLikVal )
   }

   result <- maxLik( logLik = logLik, start = start, method = method, 
      finalHessian = finalHessian, 
      yMat = yMat, xMat = xMat, randomSeed = random.seed, 
      llOneSidedGrad = oneSidedGrad, llEps = eps,
      nCoef = nCoef, nDep = nDep, ... )

   class( result ) <- c( "mvProbit", class( result ) )

   return( result )
}
