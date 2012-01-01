summary.intReg <- function(object, ...) {
   estimate <- coefTable(coef(object), stdEr(object), object$param$df)
   library(maxLik)
   s <- maxLik:::summary.maxLik(object)
   s$estimate <- estimate
   s$param <- object$param
                           # supply the additional parameters
   class(s) <- c("summary.intReg", class(s))
   return(s)
}
   
