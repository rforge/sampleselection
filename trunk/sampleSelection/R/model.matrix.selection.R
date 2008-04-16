model.matrix.selection <- function( object, part = "outcome", ... ) {

   if( !( part %in% c( "outcome", "selection" ) ) ) {
      stop( "argument 'part' must be either 'outcome' or 'selection'" )
   }

   # 2-step estimation
   if( object$method == "2step" ) {
      if( part == "selection" ) {
         result <- model.matrix( object$probit, ... )
      } else if( part == "outcome" ) {
         response <- model.frame( object$probit )[ , 1 ]
         nObs <- length( response )
         obsNames <- row.names( model.frame( object$probit ) )
         if( object$tobitType == 2 ) {
            mm <- model.matrix( object$lm, ... )
            result <- matrix( NA, nrow = nObs, ncol = ncol( mm ) )
            result[ response == 1, ] <- mm
            attributes( result )$assign <- attributes( mm )$assign
            attributes( result )$contrasts <- attributes( mm )$contrasts
            rownames( result ) <- obsNames
            colnames( result ) <- colnames( mm )
         } else if( object$tobitType == 5 ) {
            result <- list()
            mm <- list()
            mm[[ 1 ]] <- model.matrix( object$lm1, ... )
            mm[[ 2 ]] <- model.matrix( object$lm2, ... )
            for( i in 1:2 ) {
               result[[ i ]] <- matrix( NA, nrow = nObs, ncol = ncol( mm[[ i ]] ) )
               result[[ i ]][ response == ( i - 1 ), ] <- mm[[ i ]]
               attributes( result[[ i ]] )$assign <-
                  attributes( mm[[ i ]] )$assign
               attributes( result[[ i ]] )$contrasts <-
                  attributes( mm[[ i ]] )$contrasts
               rownames( result[[ i ]] ) <- obsNames
               colnames( result[[ i ]] ) <- colnames( mm[[ i ]] )
            }
         } else {
            stop( "unknown tobit type '",  object$tobitType,
               "' in object$tobitType" )
         }
      } else {
         stop( "argument 'part' must be either 'outcome' or 'selection'" )
      }
   # maximum likelihood estimation
   } else if( object$method == "ml" ) {
      stop( "the 'model.matrix' method has not yet been implemented for objects",
         " estimated by Maximum Likelihood" )
   }

   return( result )
}
