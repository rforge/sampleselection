margEff.mvProbit <- function( object, data = eval( object$call$data ),
   othDepOne = FALSE, calcVCov = FALSE,... ) {

   formula <- eval( object$call$formula )
   if( othDepOne ) {
      formula <- formula[ - 2 ]
   }

   if( calcVCov ) {
      vcov <- vcov( object )
   } else {
      vcov <- NULL
   }

   result <- mvProbitMargEff( formula = formula, 
      coef = coef( object ), vcov = vcov, data = data, ... )

   return( result )
}

