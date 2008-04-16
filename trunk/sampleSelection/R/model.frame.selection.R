model.frame.selection <- function( formula, ... ) {

   # 2-step estimation
   if( formula$method == "2step" ) {
      result <- model.frame( formula$probit, ... )
      response <- result[ , 1 ]
      attributes( result )$terms <- NULL
      obsNames <- rownames( result )
      result <- as.data.frame( cbind( result$YS, result$XS ) )
      names( result )[ 1 ] <- all.vars( formula$call$selection )[ 1 ]
      rownames( result ) <- obsNames
      if( formula$tobitType == 2 ) {
         mf <- model.frame( formula$lm, ... )
         varName <- all.vars( formula$call$outcome )[ 1 ]
         if( !varName %in% names( result ) || is.null( names( result ) ) ) {
            result[[ varName ]] <- NA
            result[[ varName ]][ response == 1 ] <- mf$YO
         }
         for( i in 1:ncol( mf$XO ) ) {
            varName <- colnames( mf$XO )[ i ]
            if( !varName %in% names( result ) ) {
               result[[ varName ]] <- NA
               result[[ varName ]][ response == 1 ] <- mf$XO[ , i ]
            }
         }
         result[[ "invMillsRatio" ]] <- NA
         result[[ "invMillsRatio" ]][ response == 1 ] <- mf[[ "imrData$IMR1" ]]
      } else if( formula$tobitType == 5 ) {
         mf <- list()
         mf[[ 1 ]] <- model.frame( formula$lm1, ... )
         mf[[ 2 ]] <- model.frame( formula$lm2, ... )
         for( eq in 1:2 ) {
            varName <- all.vars( formula$call$outcome[[ eq + 1 ]] )[ 1 ]
            if( !varName %in% names( result ) ) {
               result[[ varName ]] <- NA
            }
            result[[ varName ]][ response == ( eq - 1 ) ] <- mf[[ eq ]]$YO
            for( i in 1:( ncol( mf[[ eq ]]$XO ) - 1 ) ) {
               varName <- colnames( mf[[ eq ]]$XO )[ i ]
               if( !varName %in% names( result ) ) {
                  result[[ varName ]] <- NA
               }
               result[[ varName ]][ response == ( eq - 1 ) ] <-
                  mf[[ eq ]]$XO[ , i ]
            }
            varName <- paste( "invMillsRatio", eq, sep = "" )
            result[[ varName ]] <- NA
            result[[ varName ]][ response == ( eq - 1 ) ] <-
               mf[[ eq ]]$XO[ , "invMillsRatio" ]
         }
      } else {
         stop( "unknown tobit type '",  formula$tobitType,
            "' in formula$tobitType" )
      }
   # maximum likelihood estimation
   } else if( formula$method == "ml" ) {
      stop( "the 'model.frame' method has not yet been implemented",
         " for sample selection models estimated by Maximum Likelihood" )
   }

   return( result )
}
