margEff.mvProbit <- function( object, data = eval( object$call$data ),
   othDepOne = FALSE, ... ) {

   nCoef <- object$nDep * object$nReg

   sigma <- diag( object$nDep )
   sigma[ upper.tri( sigma ) ] <- coef( object )[ -(1:nCoef) ]
   sigma[ lower.tri( sigma ) ] <- t( sigma )[ lower.tri( sigma ) ]

   formula <- eval( object$call$formula )
   if( othDepOne ) {
      formula <- formula[ - 2 ]
   }

   result <- mvProbitMargEff( formula = formula, 
      coef = coef( object )[ 1:nCoef ], sigma = sigma, data = data, ... )

   return( result )
}

