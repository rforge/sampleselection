predict.intReg <- function ( object, newdata = NULL, type = "link", ... ) {
   
   if( is.null( newdata ) ) {
      mm <- model.matrix( object )
   }
   else {
                           # modified copy from predict.lm()
      tt <- terms( object )
      Terms <- delete.response( tt )
      m <- model.frame( Terms, newdata, xlev = object$xlevels )
      mm <- model.matrix( Terms, m )
   }
   coefnames <- colnames(mm)
   beta <- coef(object)[coefnames]
   link <- drop( mm %*% beta )
   if(type == "link") {
      ## Just the link function
      return(link)
   }
   if(type == "response") {
      ## Probabilities by intervals
      ## object$method  tells the probability distribution used
      stop("Response prediction not implemented.  What should it be?")
   }
   if(type == "linkConditional") {
      ## Expected link value given the actual interval
      if(disturbances(object) != "probit") {
         stop("conditional link prediction for ", disturbances(object),
              " disturbances not implemented")
      }
      sigma <- coef(object)["sigma"]
      if(is.null(newdata)) {
         int <- intervals(object)
         io <- intervalObs(object)
      }
      else {
         int <- intervals(newdata)
         io <- int[,"UB"] - int[,"LB"] > minIntervalWidth(object)
      }
      condE <- numeric(nrow(int))
      if(any(io)) {
         ## interval observations
         LB <- int[io,"LB"]
         UB <- int[io,"UB"]
         condE[io] <- link[io] +
             sigma*(dnorm((LB - link[io])/sigma) - dnorm((UB - link[io])/sigma))/
                 (pnorm((UB - link[io])/sigma) - pnorm((LB - link[io])/sigma))
      }
      if(!all(io)) {
         ## point observations: the conditional expectation is the mean of interval
         condE[!io] <- (int[!io,"LB"] + int[!io,"UB"])/2
      }
      return(condE)
   }
   ##   
   stop( "argument 'type' must be either 'link', 'response', or 'linkConditional'" )
}
