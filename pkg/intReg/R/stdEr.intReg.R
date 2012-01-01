stdEr.intReg <- function( x, boundaries=FALSE, ... ) {
   ## coef method.  By default, ignore the fixed boundaries.
   stde <- sqrt(diag(vcov(x, boundaries=boundaries)))
   class( stde ) <- c( "stdEr.intReg", class(stde) )
   return( stde )
}
