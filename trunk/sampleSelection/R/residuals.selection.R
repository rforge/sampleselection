residuals.selection <- function( object, part = "outcome", ... ) {

   if( !( part %in% c( "outcome", "selection" ) ) ) {
      stop( "argument 'part' must be either 'outcome' or 'selection'" )
   }

   # 2-step estimation
   if( object$method == "2step" ) {
      if( part == "selection" ) {
         result <- residuals( object$probit, ... )
      } else if( part == "outcome" ) {
         response <- model.frame( object$probit )[ , 1 ]
         result <- rep( NA, length( response ) )
         if( object$tobitType == 2 ) {
            result[ response == 1 ] <- residuals( object$lm, ... )
         } else if( object$tobitType == 5 ) {
            result[ response == 0 ] <- residuals( object$lm1, ... )
            result[ response == 1 ] <- residuals( object$lm2, ... )
         } else {
            stop( "unknown tobit type '",  object$tobitType,
               "' in object$tobitType" )
         }
         names( result ) <- row.names( model.frame( object$probit ) )
      } else {
         stop( "argument 'part' must be either 'outcome' or 'selection'" )
      }
   # maximum likelihood estimation
   } else if( object$method == "ml" ) {
      stop( "the 'residuals' method has not yet been implemented for objects",
         " estimated by Maximum Likelihood" )
   }

   return( result )
}
