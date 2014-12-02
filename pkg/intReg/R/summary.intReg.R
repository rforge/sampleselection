summary.intReg <- function(object, ...) {
   s <- NextMethod("summary", object, ...)
   estimate <- coefTable(coef(object), stdEr(object), object$param$df)
   s$estimate <- estimate
   s$param <- object$param
                           # supply the additional parameters
   class(s) <- c("summary.intReg", class(s))
   return(s)
}
   
