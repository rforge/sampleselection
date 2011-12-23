summary.intreg <- function(object, ...) {
   estimate <- coefTable(coef(object), stdEr(object), object$param$df)
   library(maxLik)
   s <- maxLik:::summary.maxLik(object)
   s$estimate <- estimate
   s$param <- object$param
                           # supply the additional parameters
   class(s) <- c("summary.intreg", class(s))
   return(s)
}

print.summary.intreg <- function(x,
                                 digits=max(3, getOption("digits") - 3),
                                 boundaries=FALSE,
                                 ...) {
   cat("--------------------------------------------\n")
   cat("Interval regression\n")
   cat( "Maximum Likelihood estimation\n" )
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(x$estimate)) {
         cat("Log-Likelihood:", logLik(x), "\n")
      }
   iPrint <- rep(TRUE, nrow(x$estimate))
   if(!boundaries) {
      iPrint[x$param$index$boundary] <- FALSE
   }
   if(!is.null(x$estimate)) {
      cat( x$param$nObs, "observations, " )
      cat( x$param$nParam, "free parameters" )
      cat( " (df = ", x$param$df, ")\n", sep="")
      printCoefmat( x$estimate[iPrint,], signif.legend = TRUE, digits = digits )
   }
   cat("--------------------------------------------\n")
   invisible( x )
}
   
