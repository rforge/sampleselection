# based on the code of "fg nu" posted at
# http://stackoverflow.com/questions/14005788/predict-function-for-heckman-model
predict.selection <- function( object, newdata = NULL,
   part = ifelse( type %in% c( "unconditional", "conditional" ),
      "outcome", "selection" ),
   type = "unconditional", ... ) {

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
   if( part == "selection" ) {
      if( type == "link" ) { 
         pred <- mXSelection %*% vBetaS
      } else if( type == "response" ) {
         pred <- pnorm( mXSelection %*% vBetaS )
      } else {
         stop( "if argument 'part' is equal to 'selection',",
            " argument 'type' must be either 'link' or 'response'" )
      }
   } else if( part == "outcome" ) {
      if( type == "unconditional" ) {
         pred <- mXOutcome %*% vBetaO
      } else if( type == "conditional" ) {
         pred <- mXOutcome %*% vBetaO +
            dnorm( temp <- mXSelection %*% vBetaS ) / pnorm( temp ) * dLambda
      } else {
         stop( "if argument 'part' is equal to 'outcome',",
            " argument 'type' must be either 'unconditional' or 'conditional'" )
      }
   } else {
      stop( "argument 'part' must be either 'selection' or 'outcome'" )
   }

   pred <- drop( pred )

   return( pred )
}
