pmvnormWrap <- function( lower = -Inf, upper = Inf, sigma, algorithm, 
   random.seed, ... ) {

   # check argument 'algorithm'
   algOkay <- TRUE
   if( is.function( algorithm ) ) {
      algResult <- do.call( algorithm, list() )
      if( ! class( algResult )[ 1 ] %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) {
         algOkay <- FALSE
      }
   } else if( is.character( algorithm ) ) {
      if( ! algorithm %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) {
         algOkay <- FALSE
      }
   } else if( ! class( algorithm ) %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) { 
      algOkay <- FALSE
   }
   if( !algOkay ) {
      stop( "argument 'algorithm' must be either one of the functions",
         " 'GenzBretz()', 'Miwa()', or 'TVPACK()'",
         " or one of the character strings",
         " \"GenzBretz\", \"Miwa\", or \"TVPACK\"" )
   }


   # checking argument 'random.seed'
   if( length( random.seed ) != 1 ) {
      stop( "argument 'random.seed' must be a single numerical value" )
   } else if( !is.numeric( random.seed ) ) {
      stop( "argument 'random.seed' must be numerical" )
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


   result <- pmvnorm( lower = lower, upper = upper, sigma = sigma,
      algorithm = algorithm, ... )

   return( result )
}
