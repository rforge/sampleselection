tobitBayes <- function( formula, data, nRep = 500, nBurn = 100, ... ) {

   if( ! "plm.dim" %in% class( data ) ) {
      stop( "argument 'data' must be a panel data frame created",
         " with 'plm.data'" )
   }

   # list to be returned at the end
   result <- list()

   # save the call that was used to call this function
   result$call <- match.call()

   # load "plm" package for paanel data estimation
   library( "plm" )

   # load "mvtnorm" package for generating multivariate normal variables
   library( "mvtnorm" )

   # detect observations that are removed due to missing values (NAs)
   obsNames <- rownames( data )
   rownames( data ) <- c( 1:nrow( data ) )
   mf <- model.frame( formula, data = data )
   validObs <- as.integer( rownames( mf ) )
   rownames( data ) <- obsNames

   # obtain and check censored dependent variable
   yCens <- unlist( mf[[1]] )
   if( min( yCens ) < 0 ) {
      stop( "the dependent variable has at least one value that is",
         " below the left limit for censoring (zero)" )
   }
   if( sum( yCens == 0 ) == 0 ) {
      stop( "the dependent variable is not censored at zero." )
   }
   obsCens <- which( yCens == 0 )

   # obtain model matrix
   xMat <- model.matrix( formula, data = data )

   # initial plm estimation to obtain starting values (step 1)
   plmStart <- plm( formula = formula, data = data, effect = "individual",
      model = "random", ... )
   # coefficients
   beta <- coef( plmStart )
   # variance of individual effects
   s2nu <- plmStart$sigma2$id
   # variance of the idiosyncratic error
   s2eps <- plmStart$sigma2$idios
                                        # individual effects
   effInd <- yCens - drop( xMat %*% beta ) - residuals( plmStart )

   # prepare formula for estimation with simulated dependent variable
   formulaSim <- formula
   formulaSim[ 2 ] <- ( tobitBayesYSim ~ . )[ 2 ]
   if( "tobitBayesYSim" %in% all.vars( formula ) ) {
      stop( "variables must not have the name 'tobitBayesYSim'" )
   }

   # matrix for all coefficients at each replication
   result$allCoef <- matrix( NA, nrow = nRep + nBurn,
      ncol = length( beta ) + 2 )
   colnames( result$allCoef ) <- c( names( beta ), "s2eps", "s2nu" )

   # MCMC replication (including burn-in)
   for( i in 1:( nBurn + nRep ) ) {
      if( i == 1 ) {
         cat( "Burn-in replications (", nBurn, "):", sep = "" )
      } else if( i == nBurn + 1 ) {
         cat( "\nMCMC replications (", nRep, "):", sep = "" )
      }
      cat( " ", i - ifelse( i <= nBurn, 0, nBurn ) )

      # simulate "latent" values to replace "censored" values
      # of the dependent variable (step 2)
      ySim <- yCens
      yFit <- drop( xMat %*% beta ) + effInd
      ySim[ obsCens ] <- NA
      for( j in obsCens ) {
         nTries <- 0
         while( is.na( ySim[ j ] ) ) {
            yCand <- rnorm( 1, mean = yFit[ j ], sd = sqrt( s2eps + s2nu ) )
            nTries <- nTries + 1
            if( yCand <= 0 ) {
               ySim[ j ] <- yCand
            } else if( nTries > 100 ) {
               ySim[ j ] <- 0
               cat( "-" )
            }
         }
      }
      data$tobitBayesYSim <- NA
      data$tobitBayesYSim[ validObs ] <- ySim

      # estimation with simulated censored variables (step 3)
      plmResult <- plm( formula = formulaSim, data = data,
         effect = "individual", model = "random", ... )
      # beta
      betaEst <- coef( plmResult )
      # covariance of beta
      betaCov <- vcov( plmResult )
      # variance of the individual effects (step 6)
      s2nuEst <- plmResult$sigma2$id
      # variance of the idiosyncratic error (step 7)
      s2epsEst <- 1 / rgamma( 1, 2, plmResult$sigma2$idios )
      # individual effects
      effInd <- yCens - drop( xMat %*% betaEst ) - residuals( plmResult )

      # draw new values for beta (step 4)
      beta <- drop( rmvnorm( 1, mean = betaEst, sigma = betaCov ) )

      # draw a new variance of the individual effects (step 6)
      s2nu <- 1 / rgamma( 1, 2, s2nuEst )

      # draw a new variance of the idiosyncratic error (step 7)
      s2eps <- 1 / rgamma( 1, 2, s2epsEst )

      # store coefficients of this replication (step 8)
      result$allCoef[ i, ] <- c( beta, s2eps, s2nu )
   }
   cat( "\n" )

   result$coef <- colMeans(
      result$allCoef[ ( nBurn + 1 ):( nBurn + nRep ), ] )
   result$coefCov <- cov(
      result$allCoef[ ( nBurn + 1 ):( nBurn + nRep ), ] )

   class( result ) <- "tobitBayes"
   return( result )
}
