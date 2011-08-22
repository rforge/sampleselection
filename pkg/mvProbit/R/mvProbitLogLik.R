mvProbitLogLik <- function( formula, coef, sigma, data,
   random.seed = 123, ... ) {

   # checking argument 'random.seed'
   if( !is.numeric( random.seed ) ) {
      stop( "argument 'random.seed' must be numerical" )
   } else if( length( random.seed ) != 1 ) {
      stop( "argument 'random.seed' must be a single numerical values" )
   }

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
