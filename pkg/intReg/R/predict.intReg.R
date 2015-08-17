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
      if(disturbancies(object) != "probit") {
         stop("conditional link prediction for ", disturbancies(object),
              " disturbancies not implemented")
      }
      sigma <- coef(object)["sigma"]
      int <- intervals(object)
      lb <- int[,"lb"]
      ub <- int[,"ub"]
      rm(int)
      condE <- link +
          sigma*(dnorm((lb - link)/sigma) - dnorm((ub - link)/sigma))/
              (pnorm((ub - link)/sigma) - pnorm((lb - link)/sigma))
      return(condE)
   }
   ##   
   stop( "argument 'type' must be either 'link', 'response', or 'linkConditional'" )
}
