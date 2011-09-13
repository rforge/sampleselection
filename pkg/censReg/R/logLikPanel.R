## log likelihood function for panel data (incl. gradients)
censRegLogLikPanel <- function( beta, yMat, xArr, left, right, nInd, nTime,
      obsBelow, obsBetween, obsAbove, nGHQ = nGHQ, ghqPoints ) {
   yMatHat <- matrix( matrix( xArr, ncol = dim( xArr )[3] ) %*%
      beta[ 1:( length( beta ) - 2 ) ], nrow = nInd, ncol = nTime )
   sigmaMu <- exp( beta[ length( beta ) - 1 ] )
   sigmaNu <- exp( beta[ length( beta ) ] )
   logLikIndMat <- matrix( NA, nrow = nInd, ncol = nGHQ )
   gradInd <- matrix( 0, nrow = nInd, ncol = length( beta ) )
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
      likGhq <- exp( logLikGhq )
      likGhqProd <- exp( logLikGhqSum )
      # gradients
      gradPartGhq <- matrix( 0, nrow = nInd, ncol = nTime )
      gradPartGhq[ obsBelow ] <-
         - dnorm( likGhqInner[ obsBelow ] ) / sigmaNu
      gradPartGhq[ obsAbove ] <-
         dnorm( likGhqInner[ obsAbove ] ) / sigmaNu
      gradPartGhq[ obsBetween ] <-
         - ddnorm( likGhqInner[ obsBetween ] ) / sigmaNu^2
      # part of gradients with respect to beta
      for( i in 1:( length( beta ) - 2 ) ) {
         gradInd[ , i ] <- gradInd[ , i ] + ghqPoints$weights[ h ] *
            likGhqProd * rowSums( gradPartGhq * xArr[ , , i ] / likGhq,
            na.rm = TRUE )
      }
      # part of gradient with respect to log( sigma_mu )
      gradInd[ , length( beta ) - 1 ] <- gradInd[ , length( beta ) - 1 ] +
         sigmaMu * ghqPoints$weights[ h ] * likGhqProd *
         rowSums( gradPartGhq * sqrt( 2 ) * ghqPoints$zeros[ h ] / likGhq )
      # part of gradient with respect to log( sigma_nu )
      gradPartGhq[ obsBelow ] <- gradPartGhq[ obsBelow ] *
         likGhqInner[ obsBelow ]
      gradPartGhq[ obsAbove ] <- - gradPartGhq[ obsAbove ] *
         likGhqInner[ obsAbove ]
      gradPartGhq[ obsBetween ] <- gradPartGhq[ obsBetween ] *
         likGhqInner[ obsBetween ] - likGhq[ obsBetween ] / sigmaNu
      gradInd[ , length( beta ) ] <- gradInd[ , length( beta ) ] +
         sigmaNu * ghqPoints$weights[ h ] * likGhqProd *
         rowSums( gradPartGhq / likGhq )
   }
   logLikInd <- rep( NA, nInd )
   for( i in 1:nInd ) {
      val <- logLikIndMat[ i, ]
      logLikInd[ i ] <- log( sum( exp( val - max( val ) ) ) ) + max( val )
   }
   ll <- logLikInd - 0.5 * log( pi )
   attr( ll, "gradient" ) <- gradInd / exp( logLikInd )
   return( ll )
}

