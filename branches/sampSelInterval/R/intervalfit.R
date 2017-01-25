intervalfit <- function(YS, XS, YO, XO, boundaries, start,
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
   ##      | 1   if  a1 < YO* <= a2 & YS = 1
   ## YO = | 2   if  a2 < YO* <= a3 & YS = 1
   ##      | ... 
   ##      \ M   if  a(M) < YO* <= a(M+1) & YS = 1
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
      sigma2 <- sqrt(exp(beta[iSigma2]))
      if( ( rho < -1) || ( rho > 1)) return(NA)
      Sigma <- matrix(c(1,-rho,-rho,1), 2, 2)
      XS.b <- drop(XS %*% betaS)
      XO.b <- drop(XO %*% betaO)
      # pre-compute the difference between the CDF of the bivariate normal
      # distribution for the interval between the boundaries
      pmvnDiff <- rep( NA, nObs )
      for( i in which(YS) ) {
         pmvnDiff[i] <-
            pmvnorm( upper = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / sigma2,
               XS.b[i] ), sigma = Sigma ) -
            pmvnorm( upper = c( ( boundaries[ YO[i] ] - XO.b[i] ) / sigma2,
               XS.b[i] ), sigma = Sigma )
      }
      loglik <- rep( NA, nObs )
      ## YS == 0, YO == NA
      loglik[YS==0] <- pnorm( -XS.b[YS==0], log.p = TRUE )
      ## YS == 1
      loglik[YS==1] <- log( pmvnDiff[YS==1] )

      ## --- gradient ---
      grad <- matrix(0, nObs, nParam)
      
      # pre-compute the difference between the PDF of the bivariate normal
      # distribution for the interval between the boundaries
      dmvnDiff <- rep( NA, nObs )
      for( i in which(YS) ) {
         dmvnDiff[i] <-
            dmvnorm( x = c( ( boundaries[ YO[i] ] - XO.b[i] ) / sigma2,
               XS.b[i] ), sigma = Sigma ) - 
            dmvnorm( x = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) / sigma2,
               XS.b[i] ), sigma = Sigma )
      }
      # gradients for the parameters for selection into policy (betaS)
      grad[YS==0, ibetaS] <-
         - dnorm( -XS.b[YS==0] ) * XS[YS==0, ] / pnorm( -XS.b[YS==0] )
      grad[YS==1, ibetaS] <- (
         pnorm( ( ( boundaries[ YO[YS==1] + 1 ] - XO.b[YS==1] ) / sigma2
            + rho * XS.b[YS==1] ) / sqrt( 1 - rho^2 ) ) -
         pnorm( ( ( boundaries[ YO[YS==1] ] - XO.b[YS==1] ) / sigma2
            + rho * XS.b[YS==1] ) / sqrt( 1 - rho^2 ) ) ) *
         dnorm( XS.b[YS==1] ) * XS[ YS==1, ] /
         pmvnDiff[YS==1]

      # gradients for the parameters for the outcome (betaO)
      grad[YS==1, ibetaO] <- (
         pnorm( ( XS.b[YS==1]
            + rho * ( ( boundaries[ YO[YS==1] + 1 ] - XO.b[YS==1] ) / sigma2 )
            ) / sqrt( 1 - rho^2 ) ) *
         dnorm( ( boundaries[ YO[YS==1] + 1 ] - XO.b[YS==1] ) / sigma2 ) - 
         pnorm( ( XS.b[YS==1]
            + rho * ( ( boundaries[ YO[YS==1] ] - XO.b[YS==1] ) / sigma2 )
            ) / sqrt( 1 - rho^2 ) ) *
         dnorm( ( boundaries[ YO[YS==1] ] - XO.b[YS==1] ) / sigma2 ) ) *
         ( -XO[ YS==1, ] / sigma2 ) /
         pmvnDiff[YS==1]

      # gradient for the correlation parameter (rho)
      grad[YS==1, iRho] <- dmvnDiff[YS==1] / pmvnDiff[YS==1] 

      # gradient for the standard deviation (sigma2)
      grad[YS==1, iSigma2] <- (
         ifelse( is.infinite( boundaries[ YO[YS==1] + 1 ] ), 0,
            pnorm( ( XS.b[YS==1] + rho *
               ( ( boundaries[ YO[YS==1] + 1 ] - XO.b[YS==1] ) / sigma2 ) ) /
                  sqrt( 1 - rho^2 ) ) *
            dnorm( ( boundaries[ YO[YS==1] + 1 ] - XO.b[YS==1] ) / sigma2 ) *
            ( ( XO.b[YS==1] - boundaries[ YO[YS==1] + 1 ] ) / sigma2^2 ) ) -
         ifelse( is.infinite( boundaries[ YO[YS==1] ] ), 0,
            pnorm( ( XS.b[YS==1] + rho *
               ( ( boundaries[ YO[YS==1] ] - XO.b[YS==1] ) / sigma2 ) ) /
                  sqrt( 1 - rho^2 ) ) *
            dnorm( ( boundaries[ YO[YS==1] ] - XO.b[YS==1] ) / sigma2 ) *
            ( ( XO.b[YS==1] - boundaries[ YO[YS==1] ] ) / sigma2^2 ) ) ) * sigma2 /
         ( pmvnDiff[YS==1] * 2 )
      
      attr(loglik, "gradient") <- grad

      return(loglik)
   } 
   

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