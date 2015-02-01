tobitTfit <- function(YS, XS, YO, XO, start,
                      print.level=0,
                      maxMethod="Newton-Raphson",
                      index=NULL,
                      binaryOutcome=FALSE,
                      ...) {
### Tobit treatment models:
### The latent variable is:
### YS* = XS'g + u
### The observables are:
###      / 1  if  YS* > 0
### YS = \ 0  if  YS* <= 0
### YO = X'b + YS bT + v
### u, v are correlated
### 
### Arguments:
### 
###  YS        binary or logical vector, 0 (FALSE) and 1 (TRUE)
###  XS              -"-                selection, should include
###              exclusion restriction
###  YO        numeric vector, outcomes
###  XO        explanatory variables for outcomes, should include YS
###  index     individual parameter indices in the parameter vector.
###            Should always be supplied but can generate here for
###            testing purposes
###  ...       additional parameters for maxLik
### 
   loglik <- function( beta) {
      betaS <- beta[iBetaS]
      betaO <- beta[iBetaO]
      sigma <- beta[iSigma]
      if(sigma <= 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
                           # check the range
      XS0.betaS <- XS0%*%betaS
                           # denoted by 'z' in the vignette
      XS1.betaS <- XS1%*%betaS
      v0 <- YO0 - XO0%*%betaO
      v1 <- YO1 - XO1%*%betaO
      sqrt1r2 <- sqrt( 1 - rho^2)
      B0 <- (-XS0.betaS - rho/sigma*v0)/sqrt1r2
      B1 <- (XS1.betaS + rho/sigma*v1)/sqrt1r2
      loglik <- numeric(nObs)
      loglik[i0] <- -1/2*log( 2*pi) - log( sigma) -
          0.5*( v0/sigma)^2 + pnorm( B0, log.p=TRUE) 
      loglik[i1] <- -1/2*log( 2*pi) -log( sigma) -
          0.5*( v1/sigma)^2 + pnorm( B1, log.p=TRUE) 
      loglik
   }
   gradlik <- function(beta) {
      ## gradient is nObs x nParam matrix
      ## components of the gradient are ordered as: g b_2, s_2, r_2, b_3,
      ## s_3, r_3
      betaS <- beta[iBetaS]
      b1 <- beta[iBetaO1]
      sigma1 <- beta[iSigma1]
      if(sigma1 <= 0) return(NA)
      rho1 <- beta[iRho1]
      if( ( rho1 < -1) || ( rho1 > 1)) return(NA)
      b2 <- beta[iBetaO2]
      sigma2 <- beta[iSigma2]
      if(sigma2 <= 0) return(NA)
      rho2 <- beta[iRho2]
      if((rho2 < -1) || (rho2 > 1)) return(NA)
                           # check the range
      XS0.betaS <- XS0%*%betaS
      XS1.betaS <- XS1%*%betaS
      XO1.b <- XO1%*%b1
      XO2.b <- XO2%*%b2
      u1 <- as.vector(YO1 - XO1.b)
      u2 <- as.vector(YO2 - XO2.b)
      sqrt1r22 <- sqrt( 1 - rho1^2)
      sqrt1r32 <- sqrt( 1 - rho2^2)
      B1 <- -(XS0.betaS + rho1/sigma1*u1)/sqrt1r22
      B2 <- (XS1.betaS + rho2/sigma2*u2)/sqrt1r32
      lambda1 <- as.vector(ifelse(B1 > -30, dnorm(B1)/pnorm(B1), -B1))
      lambda2 <- as.vector(ifelse(B2 > -30, dnorm(B2)/pnorm(B2), -B2))
                           # This is a hack in order to avoid numeric problems
      ## now the gradient itself
      gradient <- matrix(0, nObs, nParam)
      gradient[YS == 0, iBetaS] <- -lambda1*XS0/sqrt1r22
      gradient[YS == 1, iBetaS] <- lambda2*XS1/sqrt1r32
      gradient[YS == 0,iBetaO1] <- (lambda1*rho1/sigma1/sqrt1r22 + u1/sigma1^2)*XO1
      gradient[YS == 1,iBetaO2] <- (-lambda2*rho2/sigma2/sqrt1r32 + u2/sigma2^2)*XO2
      gradient[YS == 0,iSigma1] <- (-1/sigma1 + u1^2/sigma1^3
               +lambda1*rho1/sigma1^2*u1/sqrt1r22)
      gradient[YS == 1,iSigma2] <- (-1/sigma2 + u2^2/sigma2^3
               -lambda2*rho2/sigma2^2*u2/sqrt1r32)
      gradient[YS == 0,iRho1] <- -lambda1*( u1/sigma1 + rho1*XS0.betaS)/sqrt1r22^3
      gradient[YS == 1,iRho2] <- lambda2*( u2/sigma2 + rho2*XS1.betaS)/sqrt1r32^3
      gradient
   }
   hesslik <- function(beta) {
      betaS <- beta[iBetaS]
      b1 <- beta[iBetaO1]
      sigma1 <- beta[iSigma1]
      if(sigma1 <= 0) {
         return( matrix( NA, nrow = nParam, ncol = nParam ) )
      }
      rho1 <- beta[iRho1]
      if( abs( rho1 ) > 1 ) {
         return( matrix( NA, nrow = nParam, ncol = nParam ) )
      }
      b2 <- beta[iBetaO2]
      sigma2 <- beta[iSigma2]
      if(sigma2 <= 0) {
         return( matrix( NA, nrow = nParam, ncol = nParam ) )
      }
      rho2 <- beta[iRho2]
      if( abs( rho2 ) > 1 ) {
         return( matrix( NA, nrow = nParam, ncol = nParam ) )
      }
      XS0.bS <- XS0%*%betaS
      XS1.bS <- XS1%*%betaS
      u1 <- YO1 - XO1%*%b1
      u2 <- YO2 - XO2%*%b2
      sqrt1r22 <- sqrt( 1 - rho1^2)
      sqrt1r32 <- sqrt( 1 - rho2^2)
      B1 <- -(XS0.bS + rho1/sigma1*u1)/sqrt1r22
      B2 <- (XS1.bS + rho2/sigma2*u2)/sqrt1r32
      lambda1 <- ifelse(B1 > -30, dnorm(B1)/pnorm(B1), -B1)
      lambda2 <- ifelse(B2 > -30, dnorm(B2)/pnorm(B2), -B2)
                           # This is a hack in order to avoid numeric problems
      CB1 <- as.vector(ifelse(B1 > -500,
                              -exp(dnorm(B1, log = TRUE) - pnorm(B1, log.p = TRUE))*B1 -
                                  exp(2 * (dnorm(B1, log = TRUE) - pnorm(B1, log.p = TRUE))),
                              -1))
      CB2 <- as.vector(ifelse(B2 > -500,
                              -exp(dnorm(B2, log = TRUE) - pnorm(B2, log.p = TRUE))*B2 -
                                  exp(2 * (dnorm(B2, log = TRUE) - pnorm(B2, log.p = TRUE))),
                              -1))
                           # recommended by Dimitrios Rizopoulos, KULeuven
                           # This is a hack in order to avoid numerical problems.  How to do
                           # it better?  How to prove the limit value?
      l.gg <- t( XS0) %*% ( XS0 * CB1)/sqrt1r22^2 +
          t( XS1) %*% ( XS1 * CB2)/sqrt1r32^2
      l.gb1 <- -t( XS0) %*%
          ( XO1 * CB1)*rho1/sqrt1r22^2/sigma1
      l.gs2 <- -rho1/sigma1^2/sqrt1r22^2*
          t( XS0) %*% ( CB1*u1)
      l.gr2 <- t( XS0) %*%
          ( CB1*( u1/sigma1 + rho1*XS0.bS)/sqrt1r22^4 -
               lambda1*rho1/sqrt1r22^3)
      l.gb2 <- -t( XS1) %*%
          ( XO2 * CB2)*rho2/sqrt1r32^2/sigma2
      l.gs3 <- -rho2/sigma2^2/sqrt1r32^2*
          t( XS1) %*% ( CB2*u2)
      l.gr3 <- t( XS1) %*%
          ( CB2*( u2/sigma2 + rho2*XS1.bS)/sqrt1r32^4 +
               lambda2*rho2/sqrt1r32^3)
      l.b1b1 <- t( XO1) %*%
          (XO1 * ( (rho1/sqrt1r22)^2 * CB1 - 1))/sigma1^2
      l.b1s2 <- t( XO1) %*%
          ( CB1*rho1^2/sigma1^3*u1/sqrt1r22^2 -
               rho1/sigma1^2*lambda1/sqrt1r22 -
                   2*u1/sigma1^3)
      l.b1r2 <- t( XO1) %*%
          ( -CB1*( u1/sigma1 + rho1*XS0.bS)/sqrt1r22^4*rho1 +
               lambda1/sqrt1r22^3)/sigma1
      ## l.b1x3 is zero
      l.s2s2 <- sum(
          1/sigma1^2
          -3*u1*u1/sigma1^4
          + u1*u1/sigma1^4 *rho1^2/sqrt1r22^2 *CB1
          -2*lambda1* u1/sqrt1r22 *rho1/sigma1^3)
      l.s2r2 <- sum(
          ( -CB1*rho1*(u1/sigma1 + rho1*XS0.bS)/sqrt1r22 +
               lambda1)
          *u1/sigma1^2)/sqrt1r22^3
      ## l.s2x3 is zero
      l.r2r2 <- sum(
          CB1*( ( u1/sigma1 + rho1*XS0.bS)/sqrt1r22^3)^2
          -lambda1*( XS0.bS*( 1 + 2*rho1^2) + 3*rho1*u1/sigma1) /
              sqrt1r22^5
      )
      ## l.r2x3 is zero
      l.b2b2 <- t( XO2) %*%
          (XO2 * ( (rho2/sqrt1r32)^2 * CB2 - 1))/sigma2^2
      l.b2s3 <- t( XO2) %*%
          ( CB2*rho2^2/sigma2^3*u2/sqrt1r32^2 +
               rho2/sigma2^2*lambda2/sqrt1r32 - 2*u2/sigma2^3)
      l.b2r3 <- t( XO2) %*%
          ( -CB2*( u2/sigma2 + rho2*XS1.bS)/sqrt1r32^4*rho2 -
               lambda2/sqrt1r32^3)/sigma2
      l.s3s3 <- sum(
          1/sigma2^2
          -3*u2*u2/sigma2^4
          +2*lambda2* u2/sqrt1r32 *rho2/sigma2^3
          +rho2^2/sigma2^4 *u2*u2/sqrt1r32^2 *CB2)
      l.s3r3 <- -sum(
          ( CB2*rho2*(u2/sigma2 + rho2*XS1.bS)/sqrt1r32 +
               lambda2)
          *u2/sigma2^2)/sqrt1r32^3
      l.r3r3 <- sum(
          CB2*( ( u2/sigma2 + rho2*XS1.bS)/sqrt1r32^3)^2
          + lambda2*( XS1.bS*( 1 + 2*rho2^2) + 3*rho2*u2/sigma2) /
              sqrt1r32^5
      )
      hess <- array(NA, c( nParam, nParam))
      hess[iBetaS,iBetaS] <- l.gg
      hess[iBetaS,iBetaO1] <- l.gb1; hess[iBetaO1,iBetaS] <- t( l.gb1)
      hess[iBetaS,iSigma1] <- l.gs2; hess[iSigma1,iBetaS] <- t( l.gs2)
      hess[iBetaS,iRho1] <- l.gr2; hess[iRho1,iBetaS] <- t( l.gr2)
      hess[iBetaS,iBetaO2] <- l.gb2; hess[iBetaO2,iBetaS] <- t( l.gb2)
      hess[iBetaS,iSigma2] <- l.gs3; hess[iSigma2,iBetaS] <- t( l.gs3)
      hess[iBetaS,iRho2] <- l.gr3; hess[iRho2,iBetaS] <- t( l.gr3)
      hess[iBetaO1,iBetaO1] <- l.b1b1
      hess[iBetaO1,iSigma1] <- l.b1s2; hess[iSigma1,iBetaO1] <- t( l.b1s2)
      hess[iBetaO1,iRho1] <- l.b1r2; hess[iRho1,iBetaO1] <- t( l.b1r2)
      hess[iBetaO1,iBetaO2] <- 0; hess[iBetaO2,iBetaO1] <- 0
      hess[iBetaO1,iSigma2] <- 0; hess[iSigma2,iBetaO1] <- 0
      hess[iBetaO1,iRho2] <- 0; hess[iRho2,iBetaO1] <- 0
      hess[iSigma1,iSigma1] <- l.s2s2
      hess[iSigma1,iRho1] <- l.s2r2; hess[iRho1,iSigma1] <- l.s2r2
      hess[iSigma1,iBetaO2] <- 0; hess[iBetaO2,iSigma1] <- 0
      hess[iSigma1,iSigma2] <- 0; hess[iSigma2,iSigma1] <- 0
      hess[iSigma1,iRho2] <- 0; hess[iRho2,iSigma1] <- 0
      hess[iRho1,iRho1] <- l.r2r2
      hess[iRho1,iBetaO2] <- 0; hess[iBetaO2,iRho1] <- 0
      hess[iRho1,iSigma2] <- 0; hess[iSigma2,iRho1] <- 0
      hess[iRho1,iRho2] <- 0; hess[iRho2,iRho1] <- 0
      hess[iBetaO2,iBetaO2] <- l.b2b2
      hess[iBetaO2,iSigma2] <- l.b2s3; hess[iSigma2,iBetaO2] <- t( l.b2s3)
      hess[iBetaO2,iRho2] <- l.b2r3; hess[iRho2,iBetaO2] <- t( l.b2r3)
      hess[iSigma2,iSigma2] <- l.s3s3
      hess[iSigma2,iRho2] <- l.s3r3; hess[iRho2,iSigma2] <- t( l.s3r3)
      hess[iRho2,iRho2] <- l.r3r3
      hess
   }
   ## ---------------
   NXS <- ncol( XS)
   if(is.null(colnames(XS)))
      colnames(XS) <- rep("XS", NXS)
   NXO <- ncol( XO)
   if(is.null(colnames(XO)))
      colnames(XO) <- rep("XO", NXO)
   nObs <- length( YS)
   i0 <- YS==0
   i1 <- YS==1
   NO1 <- length( YS[i0])
   NO2 <- length( YS[i1])
   ## indices in for the parameter vector
   if(is.null(index)) {
      iBetaS <- 1:NXS
      iBetaO <- max(iBetaS) + seq(length=NXO)
      if(!binaryOutcome) {
         iSigma <- max(iBetaO) + 1
         iRho <- max(iSigma) + 1
      }
      else
         iRho <- max(iBetaO) + 1
      nParam <- iRho
   }
   else {
      iBetaS <- index$betaS
      iBetaO <- index$betaO
      iSigma <- index$errTerms["sigma"]
      iRho <- index$errTerms["rho"]
      nParam <- index$nParam
   }
    ## split the data by selection
    XS0 <- XS[i0,,drop=FALSE]
    XS1 <- XS[i1,,drop=FALSE]
    YO0 <- YO[i0]
    YO1 <- YO[i1]
    XO0 <- XO[i0,,drop=FALSE]
    XO1 <- XO[i1,,drop=FALSE]
    ##
    if(print.level > 0) {
        cat( "Non-participants: ", NO1,
            "; participants: ", NO2, "\n", sep="")
        cat( "Initial values:\n")
        cat("selection equation betaS:\n")
        print(start[iBetaS])
        cat("Outcome equation betaO\n")
        print(start[iBetaO])
        cat("Variance sigma\n")
        print(start[iSigma])
        cat("Correlation rho\n")
        print(start[iRho])
    }
    result <- maxLik(loglik,
#                     grad=gradlik,
#                     hess=hesslik,
                     start=start,
                     print.level=print.level,
                     method=maxMethod,
                     ...)
#   compareDerivatives(gradlik, hesslik, t0=start)
   result$tobitType <- "treatment"
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}
