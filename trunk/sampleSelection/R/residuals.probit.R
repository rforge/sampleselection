residuals.probit <- function( object, type = "response", ... ) {

   fitVal <- fitted( object )
   response <- model.frame( object )[ , 1 ]
   if( type == "response" ) {
      result <- response - fitVal
   } else if( type == "pearson" ) {
      result <- ( response - fitVal ) / sqrt( fitVal * ( 1 - fitVal ) )
   } else if( type %in% c( "deviance", "working", "partial" ) ) {
      stop( "type '", type, "' has not yet been implemented" )
   } else {
      stop( "argument 'type' must be either 'deviance', 'pearson',",
         " 'working', 'response', or 'partial'" )
   }

   result <- drop( result )
   names( result ) <- names( fitted( object ) )
   return( result )
}
