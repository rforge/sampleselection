mvProbitMargEffInternal <- function( yMat, xMat, coef, sigma,
   cond, algorithm, nGHK, eps, dummyVars,
   random.seed, ... ) {

   # number of regressors
   nReg <- ncol( xMat )

   # names of regressors
   xNames <- colnames( xMat )

   # detect dummy variables if they should be determined automatically
   if( !is.null( dummyVars ) && is.na( dummyVars[ 1 ] ) ) {
      dummyVars <- NULL
      for( i in 2:nReg ) {
         if( all( xMat[ , i ] %in% c( 0, 1, FALSE, TRUE ) ) ) {
            dummyVars <- c( dummyVars, xNames[ i ] )
         }
      }
   }

   # check dummy variables
   if( !is.null( dummyVars ) ) {
      foundDummyVars <- dummyVars %in% c( "", xNames[ - 1 ] )
      if( ! all( foundDummyVars ) ) {
         warning( "variable(s) '", 
            paste( dummyVars[ !foundDummyVars ], collapse = "', '" ),
            "' specified in argument 'dummyVars'",
            " seem(s) to be no explanatory variable(s)" )
      }
   }

   # calculate marginal effects
   for( i in 2:nReg ) {
      xMatL <- xMatU <- xMat
      isDummy <- xNames[ i ] %in% dummyVars
      if( isDummy ) {
         xMatL[ , i ] <- 0
         xMatU[ , i ] <- 1
      } else {
         xMatL[ , i ] <- xMat[ , i ] - eps / 2
         xMatU[ , i ] <- xMat[ , i ] + eps / 2
      }
      yMatL <- mvProbitExpInternal( yMat = yMat, xMat = xMatL,
         coef = coef, sigma = sigma, 
         cond = cond, algorithm = algorithm, nGHK = nGHK, 
         random.seed = random.seed, ... )
      yMatU <- mvProbitExpInternal( yMat = yMat, xMat = xMatU, 
         coef = coef, sigma = sigma, 
         cond = cond, algorithm = algorithm, nGHK = nGHK, 
         random.seed = random.seed, ... )
      margEff <- yMatU - yMatL
      if( !isDummy ) {
         margEff <- margEff / eps
      }

      # names of dependent variables
      if( i == 2 ) {
         if( !is.null( yMat ) ) {
            yNames <- colnames( yMat )
         } else {
            yNames <- paste( "y", 1:ncol( margEff ), sep = "" )
         }
      }

      # label marginal effects
      names( margEff ) <- paste( "d", yNames, "d", xNames[ i ], sep = "_" )

      if( i == 2 ) {
         result <- margEff
      } else {
         result <- cbind( result, margEff )
      }
   }

   return( result )
}
