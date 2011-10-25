margEff.mvProbit <- function( object, data = eval( object$call$data ),
   cond = FALSE, othDepOne = FALSE, dummyVars = object$dummyVars,
   calcVCov = FALSE,... ) {

   # checking argument 'data'
   if( !is.data.frame( data ) ) {
      stop( "argument 'data' must be a 'data.frame'" )
   }

   # checking argument 'cond'
   if( length( cond ) != 1 ) {
      stop( "argument 'cond' must be a single logical value" )
   } else if( !is.logical( cond ) ) {
      stop( "argument 'cond' must be a logical value" )
   }

   # checking argument 'othDepOne'
   if( length( othDepOne ) != 1 ) {
      stop( "argument 'othDepOne' must be a single logical value" )
   } else if( !is.logical( othDepOne ) ) {
      stop( "argument 'othDepOne' must be a logical value" )
   }

   # checking argument 'calcVCov'
   if( length( calcVCov ) != 1 ) {
      stop( "argument 'calcVCov' must be a single logical value" )
   } else if( !is.logical( calcVCov ) ) {
      stop( "argument 'calcVCov' must be a logical value" )
   }

   formula <- eval( object$call$formula )
   if( othDepOne || !cond ) {
      formula <- formula[ - 2 ]
   }

   if( calcVCov ) {
      vcov <- vcov( object )
   } else {
      vcov <- NULL
   }

   result <- mvProbitMargEff( formula = formula, 
      coef = coef( object ), vcov = vcov, data = data, cond = cond, 
      dummyVars = dummyVars, ... )

   return( result )
}

