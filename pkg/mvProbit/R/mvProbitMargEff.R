mvProbitMargEff <- function( formula, coef, sigma, data,
   cond = FALSE, eps = 1e-06, random.seed = 123, ... ) {

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
   if( !is.null( yMat ) ) {
      if( !is.matrix( yMat ) ) {
         stop( "either zero or at least two dependent variables",
            " must be specified in argument 'formula'",
            " (e.g. by 'cbind( y1, y2 ) ~ ...')" )
      }
   }

   # number of regressors
   nReg <- ncol( xMat )

   # names of regressors
   xNames <- colnames( xMat )

   # names of dependent variables
   if( !is.null( yMat ) ) {
      yNames <- colnames( yMat )
   } else {
      yNames <- paste( "y", 1:ncol( sigma ), sep = "" )
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

   # calculate marginal effects
   for( i in 2:nReg ) {
      xMatL <- xMatU <- xMat
      isDummy <- all( xMat[ , i ] %in% c( 0, 1, FALSE, TRUE ) )
      if( isDummy ) {
         xMatL[ , i ] <- 0
         xMatU[ , i ] <- 1
      } else {
         xMatL[ , i ] <- xMat[ , i ] - eps / 2
         xMatU[ , i ] <- xMat[ , i ] + eps / 2
      }
      yMatL <- mvProbitExp( formula, coef = coef, sigma = sigma, 
         data = as.data.frame( cbind( yMat, xMatL ) ),
         cond = cond, random.seed = random.seed, ... )
      yMatU <- mvProbitExp( formula, coef = coef, sigma = sigma, 
         data = as.data.frame( cbind( yMat, xMatU ) ),
         cond = cond, random.seed = random.seed, ... )
      margEff <- yMatU - yMatL
      if( !isDummy ) {
         margEff <- margEff / eps
      }
      names( margEff ) <- paste( "d", yNames, "d", xNames[ i ], sep = "_" )
      
      if( i == 2 ) {
         result <- margEff
      } else {
         result <- cbind( result, margEff )
      }
   }

   return( result )
}
