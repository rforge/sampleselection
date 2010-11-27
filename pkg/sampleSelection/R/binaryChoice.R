
### Binary choice models
binaryChoice <- function(formula, subset, na.action,
                   start=NULL,
                   data=sys.frame(sys.parent()),
                   x=FALSE, y=FALSE, model=FALSE,
                   method="ML",
                         cdfLower, cdfUpper, pdf, gradPdf,
                   ...) {
   ## formula: model formula, response must be either a logical or numeric vector containing only 0-s and
   ##          1-s
   ## start:      initial value of the parameters
   ## data     dataframe where the variables are defined
   ## x        whether to return model matrix
   ## y                                response vector
   ## model                            frame
   ## method   method for evaluation:
   ##          ML            maximum likelihood
   ##          model.frame   don't evaluate, only return the frame
   ## ...      further arguments for the maxLik algorithm
   ##
   ## return: a list with following components.
   ##  $results: maximisation results
   ##  $LRT:     list with two components:
   ##            LRT: likelihood ration test H0 - none of the variables significant
   ##            df:  corresponding degrees of freedom
   ##  x         model matrix (only if requested)
   ##  call      call
   ##  terms     terms
  ##  na.action  na.action used
  ##  cdfLower   lower tail of the cdf of the disturbance terms
  ##  cdfUpper   upper
  ##  pdf        pdf
  ##  gradPdf    gradient of the pdf
   loglik <- function(beta) {
      ## A generic linear index parametric binary choice log-likelihood:
      ## F: disturbance terms cdf
      ## f:                   pdf
      ## 
      ## l = sum(log(1 - F(x0'beta))) + sum(log(F(x1'beta)))
      ## dl/dbeta = sum(-f(x0'beta)/(1 - F(x0'beta)) * x0) + sum(f(x1'beta)/F(x1'beta) * x1)
      ## d2l/dbeta2 =
      ##      sum((-df(x0'beta)/dbeta*(1 - F(x0'beta)) - f^2(x0'beta))/(1 - F(x0'beta))^2 *x0'x0 +
      ##      sum((df(x1'beta)/dbeta*F(x1'beta) - f^2(x1'beta))/F(x1'beta)^2 *x1'x1
      ##
      ## We need to supply 4 specific functions:
      ## F          cdfLower    (pnorm for probit)
      ## 1 - F      cdfUpper    (pnorm(.., lower.tail=FALSE))
      ## f          pdf         (dnorm)
      ## df/dbeta   gradPdf     (-dnorm(x)*x)
      xb0 <- drop(x0 %*% beta)
      xb1 <- drop(x1 %*% beta)
      #
      F0 <- cdfUpper(xb0)
      F1 <- cdfLower(xb1)
      loglik <- numeric(length(Y))
      loglik[Y == 0] <- log(F0)
      loglik[Y == 1] <- log(F1)
      ##
      f0 <- pdf(xb0)
      f1 <- pdf(xb1)
      gradlik <- matrix(0, length(Y), length(beta))
      gradlik[Y == 0,] <- - pdf(xb0)/F0 * x0
      gradlik[Y == 1,] <- pdf(xb1)/F1 * x1
      ##
      gradf0 <- gradPdf(xb0)
      gradf1 <- gradPdf(xb1)
      hesslik <-
          t( x0) %*% ( x0 * ( -gradf0*F0 - f0*f0)/F0/F0) +
              t( x1) %*% ( x1 * ( gradf1*F1 - f1*f1)/F1/F1)
      attr(loglik, "gradient") <- gradlik
      attr(loglik, "hessian") <- hesslik
      loglik
   }
   ## Binary choice specific stuff
   cl <- match.call(call=sys.call(sys.parent()))
                           # documentation claims 'sys.call(sys.parent())' is the default
                           # value for 'call'.  But in this way it works !?
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   eval(data)
                                        # we need to eval() data here, otherwise the evaluation of the
                                        # model frame will be wrong if called from inside a function
                                        # inside a function (sorry, I can't understand it either :-(
   mf <- eval(mf, envir=parent.frame())
   if (method == "model.frame")
       return(mf)
   else if (method != "ML")
       warning("method = ", method, " is not supported. Using \"ML\"")
   mt <- attr(mf, "terms")
   Y <- model.response( mf )
   YLevels <- levels( as.factor( Y ) )
   if( length( YLevels ) != 2 ) {
      stop( "the left hand side of the 'formula' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   Y <- as.integer(Y == YLevels[ 2 ])
                                        # selection will be kept as integer internally
   X <- model.matrix(mt, mf, contrasts)
   nParam <- ncol( X)
   nObs <- length( Y)
   N1 <- sum(Y == 1)
   N0 <- nObs - N1
   if(N0 == 0 | N1 == 0) {
      stop("No variance in the response variable")
   }
   x0 <- X[Y==0,,drop=FALSE]
   x1 <- X[Y==1,,drop=FALSE]
   if(is.null(start)) {
      start <- rep( 0, nParam)
   }
   if(is.null(names(start))) {
      names(start) <- dimnames(X)[[2]]
   }
   ## Main estimation
   estimation <- maxLik(loglik, start=start,
                        method="Newton-Raphson", ...)
   ## compare.derivatives(gradlik, hesslik, t0=start)
                                        #
   ## Likelihood ratio test: H0 -- all the coefficients, except intercept
   ## are zeros.  
   ll.bar <- N0*log(N0) + N1*log(N1) - (N0 + N1)*log(N0 + N1)
                                        # log-likelihood of the H0
   LRT <- 2*(logLik(estimation) - ll.bar)
                                        # note: this assumes that the model includes constant
   result <- c(estimation,
               LRT=list(list(LRT=LRT, df=nParam-1)),
                                        # there are df-1 constraints
               param=list(list(nParam=nParam,nObs=nObs, N1=N1, N0=N0,
                               levels=YLevels)),
                           # We estimate probability for Y == levels[2]
                           # as opposite of levels[1]
               df=nObs - nParam,
               call=cl,
               terms=mt,
               x=switch(x, "1"=list(X), "0"=NULL),
               y=switch(y, "1"=list(Y), "0"=NULL),
               model=switch(model, "1"=list(mf), "0"=NULL),
               na.action=list(attr(mf, "na.action")),
               family=list(binomial(link="probit"))
                           # NA action and the removed cases
               )
   class(result) <- c("binaryChoice", class(estimation))
   result
}
