print.censReg <- function( x, logSigma = TRUE, digits = 4, ... ) {

   cat( "\n" )
   cat( "Call:\n" )
   print( x$call )
   cat( "\n" )
   cat( "Coefficients:\n" )
   print( coef( x, logSigma = logSigma ), digits = digits )
   cat( "\n" )
   invisible( x )
}
