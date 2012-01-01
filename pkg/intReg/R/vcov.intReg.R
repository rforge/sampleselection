vcov.intReg <- function( object, boundaries=FALSE, ... ) {
   ## vcov method.  By default, ignore the fixed boundaries.
   vc <- maxLik:::vcov.maxLik(object)
   if(!boundaries) {
      i <- rep(TRUE, maxLik:::nParam.maxim(object))
      i[object$param$index$boundary] <- FALSE
      vc <- vc[i,i]
   }
   class( vc) <- c( "vcov.intReg", class(vc) )
   return( vc )
}
