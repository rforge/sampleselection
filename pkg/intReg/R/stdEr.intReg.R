stdEr.intReg <- function( x, boundaries=FALSE, ...) {
   ## stdEr method.  By default, ignore the fixed boundaries.
   cat("stdEr.intReg..")
   stde <- NextMethod("stdEr", ...)
   if(!boundaries) {
      i <- rep(TRUE, nParam(x))
      i[x$param$index$boundary] <- FALSE
      stde <- stde[i]
   }
   class( stde ) <- c( "stdEr.intReg", class(stde) )
   cat("\n")
   return( stde )
}
