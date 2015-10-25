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
      }
      else {
         int <- intervals(newdata)
      }
      LB <- int[,"LB"]
      UB <- int[,"UB"]
      condE <- link +
          sigma*(dnorm((LB - link)/sigma) - dnorm((UB - link)/sigma))/
              (pnorm((UB - link)/sigma) - pnorm((LB - link)/sigma))
      return(condE)
   }
   ##   
   stop( "argument 'type' must be either 'link', 'response', or 'linkConditional'" )
}
