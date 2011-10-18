margEff.mvProbit <- function( object, data = eval( object$call$data ),
   cond = FALSE, othDepOne = FALSE, calcVCov = FALSE,... ) {

   formula <- eval( object$call$formula )
   if( othDepOne || !cond ) {
      formula <- formula[ - 2 ]
   }

   if( calcVCov ) {
      vcov <- vcov( object )
   } else {
      vcov <- NULL
   }

   result <- mvProbitMargEff( formula = formula, 
      coef = coef( object ), vcov = vcov, data = data, cond = cond, ... )

   return( result )
}

