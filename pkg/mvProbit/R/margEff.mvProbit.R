margEff.mvProbit <- function( object, data = eval( object$call$data ),
   othDepOne = FALSE, ... ) {

   formula <- eval( object$call$formula )
   if( othDepOne ) {
      formula <- formula[ - 2 ]
   }

   result <- mvProbitMargEff( formula = formula, 
      coef = coef( object ), data = data, ... )

   return( result )
}

