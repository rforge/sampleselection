intervalfit <- function(YS, XS, YO, XO, boundaries, start, AnalyticGrad,
                      weights = NULL, printLevel = 0, returnLogLikStart = FALSE,
                      maxMethod = "BHHH",
                      ...) {
   ## Fit intervaL regression model with sample selection
   ## The model is as follows:
   ## 
   ## The latent variables are:
   ## YS* = XS'BS + u1
   ## YO* = XO'BO + u2
   ## 
   ## The observables are:
   ##      / 0  if  YS* <= 0
   ## YS = \ 1  if  YS* > 0
   ## 
   ##      / NA  if  YS = 0
   ##      | 1   if  a0 < YO* <= a1 & YS = 1
   ## YO = | 2   if  a1 < YO* <= a2 & YS = 1
   ##      | ... 
   ##      \ M   if  a(M-1) < YO* <= a(M) & YS = 1
   ## 
   ##  M       number of intervals
   ##
   ##  YS      binary or logical vector, 0 (FALSE) corresponds to
   ##          YO not observed, 1 (TRUE) if observed
   ##  XS      matrix of explanatory variables for selection equation,
   ##          should include exclusion restriction
   ##  YO      outcome vector: vector of integers with values between 1 and M
   ##  XS      matrix of explanatory variables for outcomes
   ##  ...     additional parameters for maxLik
   ##  
   ## Result:
   ## Object of class 'selection', derived from 'maxLik'.
   ## Includes all the components of maxLik and additionally
   ## ...
   ## 
   loglik <- function( beta) {
      betaS <- beta[ibetaS]
      betaO <- beta[ibetaO]
      rho <- beta[iRho]
      sigma2 <- sqrt(beta[iSigma2])
      if( ( rho < -1) || ( rho > 1)) return(NA)
      Sigma <- matrix(c(1,-rho,-rho,1), 2, 2)
      XS.b <- drop(XS %*% betaS)
      XO.b <- drop(XO %*% betaO)
      loglik <- rep( NA, nObs )
      ## YS == 0, YO == NA
      loglik[YS==0] <- pnorm( -XS.b[YS==0], log.p = TRUE )
      ## YS == 1
      for( i in 1:nObs ) {
         if( YS[i] == 1 ) {
            loglik[ i ] <- log(
               pmvnorm( upper = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) /
                     sigma2, XS.b[i] ), sigma = Sigma ) -
                  pmvnorm( upper = c( ( boundaries[ YO[i] ] - XO.b[i] ) /
                        sigma2, XS.b[i] ), sigma = Sigma ) )
            # browser()
         }
      }
      ## --- gradient ---
      grad <- matrix(0, nObs, nParam)
      
      # gradients for the parameters for selection into policy (betaS)
      if(AnalyticGrad == TRUE ){
         for( i in 1:nObs ) {
            if( YS[i] == 0 ) {
               grad[i, ibetaS] <- ( dnorm( -XS.b[i] ) * XS[i] )/ pnorm( -XS.b[i] )
            } else {
               grad[i, ibetaS] <- (
                  pnorm( ( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / sigma2
                     - rho * XS.b[i] ) / sqrt( 1 - rho^2 ) ) * 
                  dnorm( XS.b[i] ) * XS[i] -
                  pnorm( ( ( boundaries[ YO[i] ] - XO.b[i] ) / sigma2
                     - rho * XS.b[i] ) / sqrt( 1 - rho^2 ) ) *
                  dnorm( XS.b[i] ) * XS[i] ) /
                  ( pmvnorm(
                     upper = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / sigma2,
                        XS.b[i] ), sigma = Sigma ) - 
                  pmvnorm(
                     upper = c( ( boundaries[ YO[i] ] - XO.b[i] ) / sigma2,
                        XS.b[i] ), sigma = Sigma ) )
            }
         }
         
         # gradients for the parameters for the outcome (betaO)
         for( i in 1:nObs ) {
            if( YS[i] == 0 ) {
               grad[i, ibetaO] <- 0
            } else {
               grad[i, ibetaO] <- (pnorm((XS.b[i] - rho * 
                     ((boundaries[ YO[i] + 1 ] - XO.b[i])/sigma2))/
                     (sqrt(1-rho^2))) * dnorm((boundaries[ YO[i] + 1 ] - 
                           XO.b[i])/sigma2) * (-XO[i]/sigma2) - 
                     pnorm((XS.b[i] - rho * ((boundaries[ YO[i] ] - XO.b[i])/
                           sigma2))/
                           (sqrt(1-rho^2))) * dnorm((boundaries[ YO[i] ] - 
                                 XO.b[i])/sigma2) * (-XO[i]/sigma2) ) /
                  (pmvnorm( upper = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / 
                        sigma2, XS.b[i] ), sigma = Sigma ) - 
                        pmvnorm( upper = c( ( boundaries[ YO[i] ] - XO.b[i] ) /
                              sigma2, XS.b[i] ), sigma = Sigma ) ) 
            }
         }
         
         # gradient for the correlation parameter (rho)
         for( i in 1:nObs ) {
            if( YS[i] == 0 ) {
               grad[i, iRho] <- 0
            } else {
               grad[i, iRho] <- (dmvnorm( x = c( ( boundaries[ YO[i] ] - 
                     XO.b[i] ) /
                     sigma2, XS.b[i] ), sigma = Sigma ) - 
                     dmvnorm( x = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) /
                           sigma2, XS.b[i] ), sigma = Sigma ) ) /
                  (pmvnorm( upper = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) /
                        sigma2, XS.b[i] ), sigma = Sigma ) - 
                        pmvnorm( upper = c( ( boundaries[ YO[i] ] - XO.b[i] ) /
                              sigma2, XS.b[i] ), sigma = Sigma ) ) 
            }
         }
         
         # gradient for the standard deviation (sigma2)
         for( i in 1:nObs ) {
            if( YS[i] == 0 ) {
               grad[i, iSigma2] <- 0
            } else {
               grad[i, iSigma2] <- (
                  pnorm( ( XS.b[i] + rho *
                     ( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / sigma2 ) ) /
                        sqrt( 1 - rho^2 ) ) *
                  dnorm( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / sigma2 ) *
                  ( ( XO.b[i] - boundaries[ YO[i] + 1 ] ) / sigma2^2 ) -
                  pnorm( ( XS.b[i] + rho *
                     ( ( boundaries[ YO[i] ] - XO.b[i] ) / sigma2 ) ) /
                        sqrt( 1 - rho^2 ) ) *
                  dnorm( ( boundaries[ YO[i] ] - XO.b[i] ) / sigma2 ) *
                  ( ( XO.b[i] - boundaries[ YO[i] ] ) / sigma2^2 ) ) /
                  ( pmvnorm( upper =
                     c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / sigma2, XS.b[i] ),
                     sigma = Sigma ) -
                  pmvnorm( upper =
                     c( ( boundaries[ YO[i] ] - XO.b[i] ) / sigma2, XS.b[i] ),
                     sigma = Sigma ) )
               # if( is.na( grad[i, iSigma2] ) ) browser()
            }
         }
         attr(loglik, "gradient") <- grad
      }
      
      return(loglik)
   } 
   
#    grad <- matrix(0, nObs, nParam)
#    r <- sqrt(1 - rho^2)
#    ## YS == 0, YO == 0
#    grad[i00,ibetaS] <- w00 * XS[i00,] * (-dnorm(-XS00.b)/lik[i00])
#    ## YS == 1, YO == 0
#    A <- dnorm(XS10.b)
#    B <- A*pnorm((XO10.b - rho*XS10.b)/r, lower.tail=FALSE)
#    grad[i10,ibetaS] <- w10 * XS[i10,]*B/lik[i10]
#    A <- dnorm(XO10.b)
#    B <- A*pnorm((XS10.b - rho*XO10.b)/r)
#    grad[i10,ibetaO] <- -w10 * XO[i10,]*B/lik[i10]
#    locmat <- -cbind(XS10.b, XO10.b)
#    pdf <- apply(locmat, 1,
#        function(x) dmvnorm(x, c(0,0), Sigma))
#    grad[i10,iRho] <- -w10 * pdf/lik[i10]
#    ## YS == 1, YO == 1
#    A <- dnorm(XS11.b)
#    B <- A*pnorm((XO11.b - rho*XS11.b)/r)
#    grad[i11,ibetaS] <- w11 * XS[i11,]*B/lik[i11]
#    A <- dnorm(XO11.b)
#    B <- A*pnorm((XS11.b - rho*XO11.b)/r)
#    grad[i11,ibetaO] <- w11 * XO[i11,]*B/lik[i11]
#    locmat <- -cbind(XS11.b, XO11.b)
#    pdf <- apply(locmat, 1,
#        function(x) dmvnorm(x, c(0,0), Sigma))
#    grad[i11,iRho] <- w11 * pdf/lik[i11]
#    ## loglik <- sum(loglik)
#    ## grad <- colSums(grad)
#    attr(loglik, "gradient") <- grad
#    return(loglik)
# }
   
   
   
   
   
   gradlik <- function(x) {
      l <- loglik(x)
      return(attr(l, "gradient"))
   }
   
   YO <- as.integer( round( YO ) )
   if( min( YO[YS==1] ) <= 0 ) {
      stop( "YO should only have strictly positive integer values" )
   }
   nInterval <- max( YO[YS==1] )
   if( length( boundaries ) != nInterval + 1 ) {
      stop( "argument 'boundaries' must have ", nInterval + 1, "elements" )
   }

   ## ---------------
   NXS <- ncol( XS )
   if(is.null(colnames(XS))) {
      colnames(XS) <- rep("XS", NXS)
   }
   NXO <- ncol( XO )
   if(is.null(colnames(XO))) {
      colnames(XO) <- rep("XO", NXO)
   }
   nObs <- length( YS )
   NO <- length( YS[YS > 0] )
   
   ## parameter indices
   ibetaS <- seq( from = 1, length.out = NXS )
   ibetaO <- seq( from = NXS+1, length.out = NXO )
   iRho <- NXS + NXO + 1
   iSigma2 <- NXS + NXO + 2
   nParam <- iSigma2
   
   # weights
   if( !is.null( weights ) ) {
      stop( "weights have not been implemented yet. Sorry!" )
   }
   
   ## output, if asked for it
   if( printLevel > 0) {
      cat("YO observed:", NO, "times; not observed:", nObs - NO,
          "times:\n")
      cat( "Number of intervals: ", nInterval, "\n" )
      print(table(YS, YO, exclude=NULL))
      cat( "Boundaries:\n")
      print(boundaries)
      cat( "Initial values:\n")
      print(start)
   }
   if( printLevel > 1) {
      cat( "Log-likelihood value at initial values:\n")
      print(loglik(start))
   }
   
   # browser()
   # # check if the likelihood values of all possible outcomes sum up to one
   # YS[] <- 0; ll0 <- loglik( start )
   # YS[] <- 1; YO[] <- 1; ll11 <- loglik( start )
   # YS[] <- 1; YO[] <- 2; ll12 <- loglik( start )
   # YS[] <- 1; YO[] <- 3; ll13 <- loglik( start )
   # all.equal( exp(ll0) + exp(ll11) + exp(ll12) + exp(ll13), rep( 1, nObs ) )
   # # check analytical derivatives
   # compareDerivatives(loglik, gradlik, t0=start )
   # range(numericGradient(loglik, t0=start)-gradlik(start))
   
   if( returnLogLikStart ) {
      return( loglik( start ) )
   }
   
   ## estimate
   result <- maxLik(loglik, 
                    start=start,
                    method=maxMethod,
                    print.level = printLevel, ... )
   result$tobitType <- "interval"
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}