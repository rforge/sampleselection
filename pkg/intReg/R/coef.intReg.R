coef.intReg <- function( object, boundaries=FALSE, ... ) {
   coefValues <- maxLik:::coef.maxLik(object)
   class( coefValues ) <- c( "coef.intReg", class(coefValues) )
   return( coefValues )
}
