residuals.probit <- function( object, type = "response", ... ) {

   if( type == "response" ) {
      result <- drop( model.frame( object )[ , 1 ] - fitted( object ) )
   } else if( type %in% c( "deviance", "pearson", "working", "partial" ) ) {
      stop( "type '", type, "' has not yet been implemented" )
   } else {
      stop( "argument 'type' must be either 'deviance', 'pearson',",
         " 'working', 'response', or 'partial'" )
   }

   names( result ) <- names( fitted( object ) )
   return( result )
}
