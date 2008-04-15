fitted.probit <- function( object, ... ) {
   mMatrix <- model.matrix( object )
   mCoef <- coef( object )
   result <- drop( pnorm( mMatrix %*% mCoef ) )
   names( result ) <- row.names( mMatrix )
   return( result )
}
