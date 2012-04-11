margEff.censReg <- function( object, ... ) {
   ## calculate marginal effects on E[y] at the mean explanatory variables
   beta <- coef( object )

   # check if the model was estimated with panel data
   isPanel <- "logSigmaMu" %in% names( beta )
   
   ## (not for panel data)
   if( !isPanel ) {
      sigma <- exp( beta[ "logSigma" ] )
      beta <- beta[ ! names(beta) %in% 
         c( "logSigma", "logSigmaMu", "logSigmaNu" ) ]
      if( length( object$xMean ) != length( beta ) ){
         warning( "cannot calculate marginal effects due to an internal error:",
                  " please contact the maintainer of this package" )
         print( beta )
         print( object$xMean )
      } else {
         xBeta <- crossprod( object$xMean, beta )
         result <- beta[ ! names( beta ) %in% c( "(Intercept)" ) ] * 
            ( pnorm( ( object$right - xBeta ) / sigma ) - 
            pnorm( ( object$left - xBeta ) / sigma ) )
         names( result ) <- 
            names( beta )[ ! names( beta ) %in% c( "(Intercept)" ) ]
      }
   } else {
      result <- NULL
   }
   return( result )
}