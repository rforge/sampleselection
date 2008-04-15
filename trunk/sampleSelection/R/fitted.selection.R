fitted.selection <- function( object, part = "outcome", ... ) {

   if( !( part %in% c( "outcome", "selection" ) ) ) {
      stop( "argument 'part' must be either 'outcome' or 'selection'" )
   }

   # 2-step estimation
   if( object$method == "2step" ) {
      if( part == "selection" ) {
         result <- fitted( object$probit, ... )
      } else if( part == "outcome" ) {
         if( object$tobitType == 2 ) {
            result <- fitted( object$lm )
         } else if( object$tobitType == 5 ) {
            response <- model.frame( object$probit )[ , 1 ]
            result <- rep( NA, length( response ) )
            result[ response == 0 ] <- fitted( object$lm1, ... )
            result[ response == 1 ] <- fitted( object$lm2, ... )
         } else {
            stop( "unknown tobit type (object$tobitType)" )
         }
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
