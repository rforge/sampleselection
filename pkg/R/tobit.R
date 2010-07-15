tobit <- function( formula, left = 0, right = Inf,
      data = sys.frame( sys.parent() ), start = NULL,
      nGHQ = 4, ... ) {

   ## checking formula
   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   } else if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   }

   ## checking limits
   # left
   if( !is.numeric( left ) ) {
      stop( "argument 'left' must be a number" )
   } else if( length( left ) != 1 ) {
      stop( "argument 'left' must be a scalar (single number)" )
   }
   # right
   if( !is.numeric( right ) ) {
      stop( "argument 'right' must be a number" )
   } else if( length( right ) != 1 ) {
      stop( "argument 'right' must be a scalar (single number)" )
   }
   # both
   if( left >= right ) {
      stop( "argument 'right' must be a larger number than argument 'left'" )
   }

   ## preparing model matrix and model response
   mc <- match.call( expand.dots = FALSE )
   m <- match( "data", names( mc ), 0 )
   mf <- mc[ c( 1, m ) ]
   mf$formula <- formula
   attributes( mf$formula ) <- NULL
   mf$na.action <- na.pass
   mf[[ 1 ]] <- as.name( "model.frame" )
   mf <- eval( mf, parent.frame() )
   mt <- attr( mf, "terms" )
   xMat <- model.matrix( mt, mf )
   xNames <- colnames( xMat )
   yVec <- model.response( mf )
   yName <- as.character( formula )[ 2 ]
   if( length( yVec ) != nrow( xMat ) ) {
      stop( "the number of observations of the endogenous variable (",
         length( yVec ), ") is not equal to the number of observations",
         " of the exogenous variables (", nrow( xMat ), ")" )
   }

   ## extract information on panel structure of data set
   isPanel <- "pdata.frame" %in% class( data )
   if( isPanel ) {
      pIndex <- attributes( data )$index
   }

   ## check if endogenous variable is within limits
   if( any( yVec < left ) ) {
      warning( "at least one value of the endogenous variable is smaller than",
         " the left limit" )
   } else if( any( yVec > right ) ) {
      warning( "at least one value of the endogenous variable is larger than",
         " the right limit" )
   }

   ## detect and remove observations with NAs, NaNs, and INFs
   validObs <- rowSums( is.na( cbind( yVec, xMat ) ) |
      is.infinite( cbind( yVec, xMat ) ) ) == 0
   yVec <- yVec[ validObs ]
   xMat <- xMat[ validObs, , drop = FALSE ]
   if( isPanel ) {
      pIndex <- pIndex[ validObs, , drop = FALSE ]
      indNames <- unique( pIndex[[ 1 ]] )  # 'names' of individuals
      nInd <- length( indNames )           # number of individuals
      timeNames <- unique( pIndex[[ 2 ]] ) # 'names' of time periods
      nTime <- length( timeNames )         # number of time periods
   }

   ## starting values
   nParam <- ncol( xMat ) + 1
   if( is.null( start ) ) {
      if( isPanel ) {
         assign( "validObs2", validObs, inherits = TRUE )
         # Random effects panel model estimation for starting values
         rEff <- plm( formula, data = data, subset = validObs2,
            effect = "individual", model = "random" )
         start <- c( coef( rEff ),
            0.5 * log( rEff$ercomp$sigma$id ),
            0.5 * log( rEff$ercomp$sigma$idios ) )
      } else {
         # OLS estimation for starting values
         ols <- lm.fit( xMat, yVec )
         start <- c( ols$coefficients,
            log( sum( ols$residuals^2 ) / length( ols$residuals ) ) )
      }
   } else {
      if( !is.numeric( start ) ) {
         stop( "argument 'start' must be numeric" )
      } else if( length( start ) != nParam ) {
         stop( "argument 'start' must have length ", nParam )
      }
   }

   if( isPanel ) {
      ## naming coefficients
      names( start ) <- c( colnames( xMat ), "logSigmaMu", "logSigmaNu" )

      ## Abscissae and weights for the Gauss-Hermite-Quadrature
      ghqPoints <- ghq( nGHQ, modified = FALSE )

      ## log likelihood function for panel data
      tobitLogLik <- function( beta ) {
         yHat <- xMat %*% beta[ 1:( length( beta ) - 2 ) ]
         sigmaMu <- exp( beta[ length( beta ) - 1 ] )
         sigmaNu <- exp( beta[ length( beta ) ] )
         ll <- rep( 0, nInd )
         for( i in 1:nInd ) {
            likInd <- 0
            for( h in 1:nGHQ ) {
               obs <- pIndex[[ 1 ]] == indNames[ i ] & yVec <= left
               tProd <- prod( 1, pnorm( ( left - yHat[ obs ] -
                  sqrt( 2 ) * sigmaMu * ghqPoints$zeros[ h ] ) / sigmaNu ) )
               obs <- pIndex[[ 1 ]] == indNames[ i ] & yVec >= right
               tProd <- prod( tProd, pnorm( ( yHat[ obs ] - right +
                  sqrt( 2 ) * sigmaMu * ghqPoints$zeros[ h ] ) / sigmaNu ) )
               obs <- pIndex[[ 1 ]] == indNames[ i ] &
                  yVec > left & yVec < right
               tProd <- prod( tProd, dnorm( ( yVec[ obs ] - yHat[ obs ] -
                  sqrt( 2 ) * sigmaMu * ghqPoints$zeros[ h ] ) / sigmaNu ) /
                  sigmaNu )
               likInd <- likInd + ghqPoints$weights[ h ] * tProd
            }
            ll[ i ] <- log( likInd / sqrt( pi ) )
         }
         return( ll )
      }
   } else {
      ## naming coefficients
      names( start ) <- c( colnames( xMat ), "logSigma" )

      ## log likelihood function for cross-sectional data
      tobitLogLik <- function( beta ) {
         yHat <- xMat %*% beta[ - length( beta ) ]
         sigma <- exp( beta[ length( beta ) ] )
         ll <- rep( NA, length( yVec ) )
         ll[ yVec <= left ] <-
            pnorm( ( left - yHat[ yVec <= left ] ) / sigma, log.p = TRUE )
         ll[ yVec > left & yVec < right ] <-
            dnorm( ( yVec - yHat )[ yVec > left & yVec < right ] / sigma,
               log = TRUE ) - log( sigma )
         ll[ yVec >= right ] <-
            pnorm( ( yHat[ yVec >= right ] - right ) / sigma, log.p = TRUE )

         ## gradients of log likelihood function for cross-sectional data
         grad <- matrix( NA, nrow = length( yVec ), ncol = length( beta ) )
         grad[ yVec <= left, ] <-
            dnorm( ( left - yHat[ yVec <= left ] ) / sigma ) /
            pnorm( ( left - yHat[ yVec <= left ] ) / sigma ) *
            cbind( - xMat[ yVec <= left, , drop = FALSE ] / sigma,
               - ( left - yHat[ yVec <= left ] ) / sigma )
         grad[ yVec > left & yVec < right, ] <-
            cbind( ( ( yVec - yHat )[ yVec > left & yVec < right ] / sigma ) *
               xMat[ yVec > left & yVec < right, , drop = FALSE ] / sigma,
               ( ( yVec - yHat )[ yVec > left & yVec < right ] / sigma )^2 - 1 )
         grad[ yVec >= right, ] <-
            dnorm( ( yHat[ yVec >= right ] - right ) / sigma ) /
            pnorm( ( yHat[ yVec >= right ] - right ) / sigma ) *
            cbind( xMat[ yVec >= right, , drop = FALSE ] / sigma,
               - ( yHat[ yVec >= right ] - right ) / sigma )
         attr( ll, "gradient" ) <- grad
         return( ll )
      }
   }

   result <- maxLik( tobitLogLik, start = start, ... )

   class( result ) <- c( "tobit", class( result ) )
   return( result )
}

