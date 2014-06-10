# based on the code of "fg nu" posted at
# http://stackoverflow.com/questions/14005788/predict-function-for-heckman-model
predict.selection <- function( object, newdata = NULL, type = "uncond", ... ) {

   if( is.null( newdata ) ) {
      # regressor matrix for the selection equation
      mXSelection <- model.matrix( object, part = "selection" )
      
      # regressor matrix for the outcome equation
      mXOutcome = model.matrix( object, part = "outcome" )
   } else {
      # construct the Formula object
      tempS <- evalq( object$call$selection )
      tempO <- evalq( object$call$outcome )
      
      FormHeck <- as.Formula(
         paste0( tempO[2], '|', tempS[2], '~', tempO[3], '|', tempS[3] ) )
      
      # regressor matrix for the selection equation
      mXSelection <- model.matrix( FormHeck, data = newdata, rhs = 2 )
      
      # regressor matrix for the outcome equation
      mXOutcome = model.matrix(FormHeck, data = newdata, rhs = 1 )
   }

   # indices of the various parameters in selectionObject$estimate
   vIndexBetaS <- object$param$index$betaS
   vIndexBetaO <- object$param$index$betaO
   vIndexErr <- object$param$index$errTerms
   
   # get the estimates
   vBetaS <- object$estimate[ vIndexBetaS ]
   vBetaO <- object$estimate[ vIndexBetaO ]
   
   dLambda <- object$estimate[ vIndexErr['rho'] ] *
      object$estimate[ vIndexErr['sigma'] ]
   
   # depending on the type of prediction requested, return
   # TODO allow the return of multiple prediction types
   if( type == "link" ) { 
      pred <- mXSelection %*% vBetaS
   } else if( type == "prob" ) {
      pred <- pnorm( mXSelection %*% vBetaS )
   } else if( type == "uncond" ) {
      pred <- mXOutcome %*% vBetaO
   } else if( type == "cond" ) {
      pred <- mXOutcome %*% vBetaO +
         dnorm( temp <- mXSelection %*% vBetaS ) / pnorm( temp ) * dLambda
   } else {
      stop( "argument 'type' must be either 'link', 'prob', 'uncond', or 'cond'" )
   }

   pred <- drop( pred )

   return( pred )
}
