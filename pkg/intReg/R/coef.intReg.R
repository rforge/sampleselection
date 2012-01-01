coef.intReg <- function( object, boundaries=FALSE, ... ) {
   ## coef method.  By default, ignore the fixed boundaries.
   coefValues <- maxLik:::coef.maxLik(object)
   if(!boundaries) {
      i <- rep(TRUE, length(coefValues))
      i[object$param$index$boundary] <- FALSE
      coefValues <- coefValues[i]
   }
   class( coefValues ) <- c( "coef.intReg", class(coefValues) )
   return( coefValues )
}
