pmvnormWrap <- function( lower = -Inf, upper = Inf, sigma, algorithm, ... ) {

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

   result <- pmvnorm( lower = lower, upper = upper, sigma = sigma,
      algorithm = algorithm, ... )

   return( result )
}
