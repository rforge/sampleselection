## log likelihood function for panel data (incl. gradients)
censRegLogLikPanel <- function( beta, yMat, xMat, xArr, left, right, nInd, nTime,
      obsBelow, obsBetween, obsAbove, nGHQ = nGHQ, ghqPoints ) {
   yMatHat <- matrix( matrix( xArr, ncol = ncol( xMat ) ) %*%
      beta[ 1:( length( beta ) - 2 ) ], nrow = nInd, ncol = nTime )
   sigmaMu <- exp( beta[ length( beta ) - 1 ] )
   sigmaNu <- exp( beta[ length( beta ) ] )
   likInd <- rep( 0, nInd )
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
      likGhq <- matrix( 1, nrow = nInd, ncol = nTime )
      likGhq[ obsBelow | obsAbove ] <-
         pnorm( likGhqInner[ obsBelow | obsAbove ] )
      likGhq[ obsBetween ] <-
         dnorm( likGhqInner[ obsBetween ] ) / sigmaNu
      likGhqProd <- apply( likGhq, 1, prod )
      likInd <- likInd + ghqPoints$weights[ h ] * likGhqProd
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
   ll <- log( likInd / sqrt( pi ) )
   attr( ll, "gradient" ) <- gradInd / likInd
   return( ll )
}

