intervalfit <- function(YS, XS, YO, XO, boundaries, start,
                      weights = NULL, printLevel = 0,
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
   ##      | 0   if  a0 < YO* <= a1 & YS = 1
   ## YO = | 1   if  a1 < YO* <= a2 & YS = 1
   ##      | ... 
   ##      \ M-1 if  a(M-1) < YO* <= a(M) & YS = 1
   ## 
   ##  M       number of intervals
   ##
   ##  YS      binary or logical vector, 0 (FALSE) corresponds to
   ##          YO not observed, 1 (TRUE) if observed
   ##  XS      matrix of explanatory variables for selection equation,
   ##          should include exclusion restriction
   ##  YO      outcome vector: vector of integers with values between 0 and M-1
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
               pmvnorm( upper = c( ( boundaries[ YO[i] + 2 ] - XO.b[i] ) /
                     sigma2, XS.b[i] ), sigma = Sigma ) -
               pmvnorm( upper = c( ( boundaries[ YO[i] + 1 ] - XO.b[i] ) /
                     sigma2, XS.b[i] ), sigma = Sigma ) )
            # browser()
         }
      }
      ## --- gradient ---
      grad <- matrix(0, nObs, nParam)
      # attr(loglik, "gradient") <- grad
      return(loglik)
   }
   gradlik <- function(x) {
      l <- loglik(x)
      return(attr(l, "gradient"))
   }
   
   YO <- as.integer( round( YO ) )
   if( min( YO[YS==1] ) < 0 ) {
      stop( "YO should only have positive values" )
   }
   nInterval <- max( YO[YS==1] ) + 1
   if( length( boundaries ) != nInterval + 1 ) {
      stop( "agrument 'boundaries' must have ", nInterval + 1, "elements" )
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
   # compareDerivatives(loglik, gradlik, t0=start )
   # range(numericGradient(loglik, t0=start)-gradlik(start))
   
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
