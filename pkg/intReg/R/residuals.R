residuals.intReg <- function(object, ... ) {
   ## compute the generalized residuals:
   ## E(eps|x, interval)
   ##
   if(disturbances(object) != "probit") {
      stop("residuals for ", disturbances(object),
           " disturbances not implemented")
   }
   int <- intervals(object)
   io <- intervalObs(object)
   l <- int[,"LB"]
   u <- int[,"UB"]
   beta <- coef(object)[object$param$index$beta]
                           # remove 'sigma' and other nuisance params
   sigma <- coef(object)["sigma"]
   mm <- model.matrix( object )
   link <- drop( mm %*% beta )
   eps <- numeric(nrow(int))
   if(any(io)) {
      ## interval observations
      LB <- int[io,"LB"]
      UB <- int[io,"UB"]
      eps[io] <- sigma*(dnorm((LB - link[io])/sigma) - dnorm((UB - link[io])/sigma))/
          (pnorm((UB - link[io])/sigma) - pnorm((LB - link[io])/sigma))
      }
   if(!all(io)) {
      ## point observations: the conditional expectation is the mean of interval
      eps[!io] <- (int[!io,"LB"] + int[!io,"UB"])/2 - link[!io]
   }
   return(eps)
}
