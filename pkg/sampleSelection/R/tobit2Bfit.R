tobit2Bfit <- function(YS, XS, YO, XO, start,
                      print.level=0,
                      maxMethod="BHHH",
                      ...) {
   ## Fit 2-dimensional sample selection models where the outcome
   ## variable is binary.
   ## The model is as follows (following Amemiya 1985):
   ## 
   ## The latent variables are:
   ## YS* = XS'g + u1
   ## YO* = X'b2 + u2
   ## 
   ## The observables are:
   ##      / 1  if  YS* > 0
   ## YS = \ 0  if  YS* <= 0
   ##      / 1(YO* > 0)  if  YS = 0
   ## YO = \ 0           if  YS = 1
   ## 
   ##  YS      binary or logical vector, 0 (FALSE) corresponds to
   ##          YO not observed, 1 (TRUE) if observed
   ##  XS      matrix of explanatory variables for selection equation,
   ##          should include exclusion restriction
   ##  YO      binary or logical outcome vector
   ##  XS      matrix of explanatory variables for outcomes
   ##  ...     additional parameters for maxLik
   ##  
   ## Result:
   ## Object of class 'tobit2', derived from 'maxLik'.
   ## Includes all the components of maxLik and additionally
   ## twoStep   Results for Heckman two-step estimator, used for initial values
   ## 
   loglik <- function( beta) {
      g <- beta[ibetaS]
      b <- beta[ibetaO]
      rho <- beta[iRho]
      Sigma <- matrix(c(1,rho,rho,1), 2, 2)
      if( ( rho < -1) || ( rho > 1)) return(NA)
      XS00.g <- XS[i00,] %*% g
      XS10.g <- XS[i10,] %*% g
      XS11.g <- XS[i11,] %*% g
      XO10.b <- XO[i10,] %*% b
      XO11.b <- XO[i11,] %*% b
      loglik <- numeric(nObs)
      loglik[i00] <- pnorm(-XS00.g, log=TRUE)
      uppermat <- -cbind(XS10.g, XO10.b)
      f2 <- apply(uppermat, 1,
                  function(x) pmvnorm(lower=x, corr=Sigma))
      loglik[i10] <- log(pnorm(-XS10.g, lower.tail=FALSE) - f2)
      uppermat <- -cbind(XS11.g, XO11.b)
      f2 <- apply(uppermat, 1,
                  function(x) pmvnorm(lower=x, corr=Sigma))
      loglik[i11] <- log(f2)
      loglik
   }
    ## ---------------
    NXS <- ncol( XS)
    if(is.null(colnames(XS)))
        colnames(XS) <- rep("XS", NXS)
    NXO <- ncol( XO)
    if(is.null(colnames(XO)))
        colnames(XO) <- rep("XO", NXO)
   nObs <- length( YS)
    NO <- length( YS[YS > 0])
   ## selection indices
   i00 <- YS == 0
   i10 <- (YS == 1) & (YO == 0)
   i11 <- (YS == 1) & (YO == 1)
   ## parameter indices
   ibetaS <- 1:NXS
   ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
   iRho <- tail(ibetaO, 1) + 1
   ## output, if asked for it
   if( print.level > 0) {
      cat("YO observed:", NO, "times; not observed:", nObs - NO,
          "times:\n")
      print(table(YS, YO, exclude=NULL))
      cat( "Initial values:\n")
      print(start)
   }
   ## pre-calculate a few values:
   XS0 <- XS[YS==0,,drop=FALSE]
   XS1 <- XS[YS==1,,drop=FALSE]
   YO[is.na(YO)] <- 0
   YO1 <- YO[YS==1]
   XO1 <- XO[YS==1,,drop=FALSE]
   N0 <- sum(YS==0)
   N1 <- sum(YS==1)
   ## estimate
##    compareDerivatives(loglik, gradlik, t0=start)
##    stop()
   library(mvtnorm)
   result <- maxLik(loglik,
                    start=start,
                    method=maxMethod,
                    print.level=print.level, ...)
   result$tobitType <- 2
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}
