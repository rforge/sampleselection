## log likelihood function for panel data (incl. gradients)
censRegLogLikPanel <- function( beta, yMat, xArr, left, right, nInd, nTime,
      obsBelow, obsBetween, obsAbove, nGHQ = nGHQ, ghqPoints ) {
   yMatHat <- matrix( matrix( xArr, ncol = dim( xArr )[3] ) %*%
      beta[ 1:( length( beta ) - 2 ) ], nrow = nInd, ncol = nTime )
   sigmaMu <- exp( beta[ length( beta ) - 1 ] )
   sigmaNu <- exp( beta[ length( beta ) ] )
   logLikIndMat <- matrix( NA, nrow = nInd, ncol = nGHQ )
   grad1LogIndArr <- array( NA, c( nInd, length( beta ), nGHQ ) )
   grad2IndArr <- array( NA, c( nInd, length( beta ), nGHQ ) )
   for( h in 1:nGHQ ) {
      likGhqInner <- matrix( NA, nrow = nInd, ncol = nTime )
      likGhqInner[ obsBelow ] <-
         ( left - yMatHat[ obsBelow ] - sqrt( 2 ) * sigmaMu *
            ghqPoints$zeros[ h ] ) / sigmaNu
      likGhqInner[ obsAbove ] <-
         ( yMatHat[ obsAbove ] - right + sqrt( 2 ) * sigmaMu *
            ghqPoints$zeros[ h ] ) / sigmaNu
      likGhqInner[ obsBetween ] <-
         ( yMat[ obsBetween ] - yMatHat[ obsBetween ] -
            sqrt( 2 ) * sigmaMu * ghqPoints$zeros[ h ] ) / sigmaNu
      logLikGhq <- matrix( 0, nrow = nInd, ncol = nTime )
      logLikGhq[ obsBelow | obsAbove ] <-
         pnorm( likGhqInner[ obsBelow | obsAbove ], log.p = TRUE )
      logLikGhq[ obsBetween ] <-
         dnorm( likGhqInner[ obsBetween ], log = TRUE ) - log( sigmaNu )
      logLikGhqSum <- apply( logLikGhq, 1, sum )
      logLikIndMat[ , h ] <- log( ghqPoints$weights[ h ] ) + logLikGhqSum
      # gradients
      gradPartGhq <- matrix( 0, nrow = nInd, ncol = nTime )
      gradPartGhq[ obsBelow ] <-
         - dnorm( likGhqInner[ obsBelow ] ) / sigmaNu
      gradPartGhq[ obsAbove ] <-
         dnorm( likGhqInner[ obsAbove ] ) / sigmaNu
      gradPartGhq[ obsBetween ] <-
         dnorm( likGhqInner[ obsBetween ] ) * likGhqInner[ obsBetween ] / sigmaNu^2
      # part of gradients with respect to beta
      for( i in 1:( length( beta ) - 2 ) ) {
         grad1LogIndArr[ , i, h ] <- log( ghqPoints$weights[ h ] ) +
            logLikGhqSum
         grad2IndArr[ , i, h ] <- 
            rowSums( gradPartGhq * xArr[ , , i ] / exp( logLikGhq ), 
            na.rm = TRUE )
      }
      # part of gradient with respect to log( sigma_mu )
      grad1LogIndArr[ , length( beta ) - 1, h ] <-
         log( sigmaMu ) + log( ghqPoints$weights[ h ] ) + logLikGhqSum
      grad2IndArr[ , length( beta ) - 1, h ] <-
         rowSums( gradPartGhq * sqrt( 2 ) * ghqPoints$zeros[ h ] / exp( logLikGhq ) )
      # part of gradient with respect to log( sigma_nu )
      gradPartGhq[ obsBelow ] <- gradPartGhq[ obsBelow ] *
         likGhqInner[ obsBelow ]
      gradPartGhq[ obsAbove ] <- - gradPartGhq[ obsAbove ] *
         likGhqInner[ obsAbove ]
      gradPartGhq[ obsBetween ] <- gradPartGhq[ obsBetween ] *
         likGhqInner[ obsBetween ] - exp( logLikGhq[ obsBetween ] ) / sigmaNu
      grad1LogIndArr[ , length( beta ), h ] <-
         log( sigmaNu ) + log( ghqPoints$weights[ h ] ) + logLikGhqSum
      grad2IndArr[ , length( beta ), h ] <-
         rowSums( gradPartGhq / exp( logLikGhq ) )
   }
   logLikInd <- rep( NA, nInd )
   gradInd <- matrix( NA, nrow = nInd, ncol = length( beta ) )
   for( i in 1:nInd ) {
      val <- logLikIndMat[ i, ]
      logLikInd[ i ] <- log( sum( exp( val - max( val ) ) ) ) + max( val )
      for( j in 1:length( beta ) ) {
         gradInd[ i, j ] <- sum( 
            exp( grad1LogIndArr[ i, j, ] - logLikInd[ i ] ) * 
            grad2IndArr[ i, j, ] )
      }
   }
   ll <- logLikInd - 0.5 * log( pi )
   attr( ll, "gradient" ) <- gradInd
   return( ll )
}

