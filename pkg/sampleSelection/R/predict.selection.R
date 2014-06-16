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
      mXOutcome <- model.matrix( object, part = "outcome" )

      # remove inverse Mills ratio
      if( object$method == "2step" ) {
         mXOutcome <- mXOutcome[ , -ncol( mXOutcome ) ]
      }
      
   } else {
      # construct the Formula object
      tempS <- eval( object$call$selection )
      tempO <- eval( object$call$outcome )
      
      formS <- as.formula( tempS )
      if( object$tobitType == 2 ) {
         formO <- as.formula( tempO )
      } else if( object$tobitType == 5 ) {
         formO <- as.formula( tempO[[1]] )
      } else {
         stop( "internal error: unknown tobitType '", object$tobitType,
            "' Please contact the maintainer of the sampleSelection package" )
      }

      # regressor matrix for the selection equation
      mfS <- model.frame( formS, data = newdata, na.action = na.pass )
      mXSelection <- model.matrix( formS, mfS )
      
      # regressor matrix for the outcome equation
      mfO <- model.frame( formO, data = newdata, na.action = na.pass )
      mXOutcome <- model.matrix( formO, mfO )
   }

   # indices of the various parameters in selectionObject$estimate
   vIndexBetaS <- object$param$index$betaS
   vIndexBetaO <- object$param$index$betaO
   
   # get the estimates
   vBetaS <- coef( object )[ vIndexBetaS ]
   vBetaO <- coef( object )[ vIndexBetaO ]
   
   dLambda <- coef( object )[ "rho" ] * coef( object )[ "sigma" ]

   # remove coefficient of inverse Mills ratio
   if( object$method == "2step" ) {
      vBetaO <- vBetaO[ names( vBetaO ) != "invMillsRatio" ]
   }
   
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
         mXOutcome <- mXOutcome[
            rownames( mXOutcome ) %in% rownames( mXSelection ), ]
         mXSelection <- mXSelection[
            rownames( mXSelection ) %in% rownames( mXOutcome ), ]
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
