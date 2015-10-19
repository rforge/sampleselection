
### Based on polr() function in MASS (originally developed for categorical
### wage data in Social Capital Benchmark Survey)
intReg <- function(formula, start, boundaries,
                   ...,
                 contrasts = NULL, Hess = FALSE,
                 model = TRUE,
                 method = c("probit", "logistic", "cloglog", "cauchit",
                 "model.frame"),
                   minIntervalWidth=10*sqrt(.Machine$double.eps),
                     print.level=0,
                   data, subset, weights, na.action,
                   iterlim=100)
{
   ## Estimates the interval regression model by Maximum Likelihood
   ## 
   ## minIntervalWidth    minimum interval width where still to consider it an interval
   ##                     if bounds are no further apart than that
   ##                     it is considered a point estimate.
   ##                
    logit <- function(p) log(p/(1 - p))
    ## log likelihood
    loglik <- function(theta) {
       beta <- theta[iBeta]
        zeta <- theta[iBoundaries]
        sigma <- theta[iSigma]
        eta <- offset
       if (nBeta > 0) {
          eta <- eta + drop(x %*% beta)
       }
       ll <- numeric(nrow(x))
       ## interval observations
       if(nIntervalObs > 0) {
          ub <- zeta[boundaryInterval[iIntervalObs] + 1]
          lb <- zeta[boundaryInterval[iIntervalObs]]
          Pr <- pfun((ub - eta[iIntervalObs])/sigma) -
              pfun((lb - eta[iIntervalObs])/sigma)
          if(any(Pr <= 0))
             return(NA)
          names(Pr) <- NULL
                           # individual probability to be in interval
                           # (zeta[y+1], zeta[y]])
          ll[iIntervalObs] <- wt[iIntervalObs] * log(Pr)
       }
       ## point observations
       if(nPointObs > 0) {
          ll[!iIntervalObs] <- wt[!iIntervalObs]*
              dfun((y[!iIntervalObs,1] - eta[!iIntervalObs])/sigma, log=TRUE) -
                  log(sigma)
                           # use lower bound for point estimate (should be equal to
                           # upper bound anyway)
       }
       ## ll <- sum(ll)
       return(ll)
     }
     ## gradient
     gradlik <- function(theta)
     {
         beta <- theta[iBeta]
         zeta <- theta[iBoundaries]
         sigma <- theta[iSigma]
         eta <- offset
         if(nBeta > 0)
             eta <- eta + drop(x %*% beta)
         ##
         grad <- matrix(0, nrow(x), length(theta))
         ## Gradient of interval observations
         if(nIntervalObs > 0) {
            normArgUB <- (zeta[boundaryInterval[iIntervalObs] +1] -
                              eta[iIntervalObs])/sigma
            normArgLB <- (zeta[boundaryInterval[iIntervalObs]] -
                              eta[iIntervalObs])/sigma
            Pr <- pfun(normArgUB) - pfun(normArgLB)
            pUB <- dfun(normArgUB)
            pLB <- dfun(normArgLB)
            dg.dbeta <-
                            # d loglik / d beta
               if(nBeta > 0)
                  x[iIntervalObs,] * (wt[iIntervalObs]*(pLB - pUB)/Pr/sigma)
               else
                  NULL
            sUB <- pUB*normArgUB
            sUB[is.infinite(normArgUB)] <- 0
                            # if a boundary is Inf, we get 0*inf type of NaN
            sLB <- pLB*normArgLB
            sLB[is.infinite(normArgLB)] <- 0
            dg.dsigma <- -(sUB - sLB)*(wt[iIntervalObs]/Pr/sigma)
            grad[iIntervalObs, iBeta] <- dg.dbeta
            grad[iIntervalObs, iSigma] <- dg.dsigma
         }
         ## Gradient of point observations
         if(nPointObs > 0) {
            normArg <- (y[!iIntervalObs,1] - eta[!iIntervalObs])/sigma
            if(method == "probit") {
               dg.dbeta <- x[!iIntervalObs,] * normArg/sigma
               dg.dsigma <- (normArg^2 - 1)/sigma
            }
            else {
               dg.dbeta <- -x[!iIntervalObs,] * d2fun(normArg)/dfun(normArg)/sigma
               dg.dsigma <- d2fun(normArg)/dfun(normArg)*normArg/sigma
            }
            i <- is.infinite(normArg)
            ##
            grad[!iIntervalObs, iBeta] <- dg.dbeta
            grad[!iIntervalObs, iSigma] <- dg.dsigma
         } 
         ## grad <- colSums(grad)
         return(grad)
      }
     ## ---------- main function -------------
     cl <- match.call(expand.dots = FALSE)
     method <- match.arg(method)
     if(is.matrix(eval.parent(cl$data)))
         cl$data <- as.data.frame(data)
     cl$start <- cl$Hess <- cl$method <- cl$model <- cl$boundaries <- cl$... <-
        cl$print.level <- cl$iterlim <- cl$minIntervalWidth <- NULL
     cl[[1]] <- as.name("model.frame")
     mf <- eval.parent(cl)
     if (method == "model.frame")
         return(mf)
     mt <- attr(mf, "terms")
     ## Select the correct model
     pfun <- switch(method, logistic = plogis, probit = pnorm,
                    cloglog = pgumbel, cauchit = pcauchy)
     dfun <- switch(method, logistic = dlogis, probit = dnorm,
                    cloglog = dgumbel, cauchit = dcauchy)
     d2fun <- switch(method,
                     probit = function(x) -x*dnorm(x),
                     function(x) stop("Point observations for ", method,
                                      " disturbances not implemented\n")
                     )
                            # derivative of the density
     ##
     x <- model.matrix(mt, mf, contrasts)
     xint <- match("(Intercept)", colnames(x), nomatch=0)
     nObs <- nrow(x)
     nBeta <- ncol(x)
     cons <- attr(x, "contrasts") # will get dropped by subsetting
     wt <- model.weights(mf)
     if(!length(wt)) {
        wt <- rep(1, nObs)
     }
     offset <- model.offset(mf)
    if(length(offset) <= 1) {
       offset <- rep(0, nObs)
    }
     y <- model.response(mf)
     if(is.matrix(y)) {
        ## Use the intervals, given for each observation.  We have interval regression and the interval boundaries are fixed.
        ## Save boundaries as a sequence of pairs of L,B boundaries
        if(ncol(y) != 2) {
           stop("response must be a factor or Nx2 matrix of boundaries")
        }
        ordered <- FALSE
        dimnames(y) <- NULL
        lowerBound <- y[,1]
        upperBound <- y[,2]
        iIntervalObs <- abs(upperBound - lowerBound) > minIntervalWidth
                           # these are interval observations
        iIntervalObs[is.infinite(lowerBound) & is.infinite(upperBound)] <- FALSE
                           # treat (Inf,Inf) intervals as point obs
        nIntervalObs <- sum(iIntervalObs)
        nPointObs <- sum(!iIntervalObs)
        ## in case of interval regression, we have to construct a set of intervals and pack them correctly to the
        ## parameter vector
        intervals <- sets::set()
        for(i in seq(length=nrow(x))) {
           intervals <- intervals | c(lowerBound[i], upperBound[i])
        }
        ## Now make the set to a list to have ordering
        intervals <- as.list(intervals)
        ## Now find which interval each observation falls into
        yInt <- boundaryInterval <- numeric(nrow(x))
                           # which interval observation falls to 
                           # yInt:                 in terms of ordered intervals
                           # boundaryInterval   in terms of ordered boundaries
        for(i in seq(along=intervals)) {
           j <- lowerBound == intervals[[i]][1] & upperBound == intervals[[i]][2]
                           # Note: y and boundaryInterval for point estimates will
                           # also be written but not used later
           boundaryInterval[j] <- 1 + 2*(i - 1)
           yInt[j] <- i
        }
        boundaries <- unlist(intervals)
                           # boundaries as a vector (note the joint boundaries are twice
        names(boundaries) <- paste(c("L", "U"), rep(seq(along=intervals), each=2))
        nInterval <- length(boundaries) - 1
     }
     else {
        ## response is given as an ordered factor, boundaries must be given separately.
        ## Save them as a vector of boundaries, all numbers (except the first, last) represent the upper boundary of
        ## the smaller interval and the lower boundary of the upper interval at the same time.
        if(!is.factor(y)) {
           stop("response must be a factor or Nx2 matrix of boundaries")
        }
        lev <- levels(y)
        if(length(lev) <= 2)
            stop("response must have 3 or more levels")
        yInt <- unclass(y)
                           # which interval observation falls to: we keep a separate
                           # variable to preserve the original
        nInterval <- length(lev)
        nIntervalObs <- 0
                           # should handle point observations here as well ....
        if(missing(boundaries))
            ordered <- TRUE
        else {
            ordered <- FALSE
            intervals <- vector("list", length(boundaries) - 1)
            for(i in seq(length=length(boundaries) - 1)) {
               intervals[[i]] <- c(boundaries[i], boundaries[i+1])
            }
            boundaryInterval <- yInt
                           # y falls inbetween boundaries 'boundaryInterval' and
                           # 'boundaryInterval + 1'
         }
        if(is.null(names(boundaries)))
            names(boundaries) <- paste("Boundary", seq(along=boundaries))
        iIntervalObs <- rep(TRUE, length(y))
        nIntervalObs <- sum(iIntervalObs)
        nPointObs <- sum(!iIntervalObs)
     }
    if(nIntervalObs > 0) {
        Y <- matrix(0, nObs, nInterval + 1)
        .polrY1 <- col(Y) == boundaryInterval + 1
        .polrY2 <- col(Y) == boundaryInterval 
                                         # .polr are markers for which interval the
                                         # boundaryInterval falls to
     }
    ## Remove unsuitable observations
    iExclude <- !complete.cases(x)
    LB <- boundaries[boundaryInterval]
    UB <- boundaries[boundaryInterval + 1]
    iExclude <- iExclude | !complete.cases(LB, UB)
    iExclude <- iExclude | (is.infinite(LB) & is.infinite(UB))
                           # exclude (Inf, Inf) intervals
#    y <- y[!iExclude,,drop=FALSE]
    boundaryInterval <- boundaryInterval[!iExclude]
    yInt <- yInt[!iExclude]
    x <- x[!iExclude,,drop=FALSE]
    wt <- wt[!iExclude]
    offset <- offset[!iExclude]
    iIntervalObsInf <- iIntervalObs
    iIntervalObs <- iIntervalObs[!iExclude]
     ## starting values
     iBeta <- seq(length=ncol(x))
                            # coefficients
     iBoundaries <- nBeta + seq(along=boundaries)
                            # boundaries
     iSigma <- max(iBeta, iBoundaries) + 1
                            # standard deviation
     if(missing(start)) {
        start <- numeric(max(iBeta, iBoundaries, iSigma))
                            # +1 for the error variance 'sigma'
        activePar <- logical(length(start))
        if(ordered) {
           ## try logistic/probit regression on 'middle' cut
           q1 <- nInterval %/% 2
           y1 <- (y > q1)
           fit <-
               switch(method,
                      "logistic"= glm.fit(x, y1, wt, family = binomial(), offset = offset),
                      "probit" = glm.fit(x, y1, wt, family = binomial("probit"), offset = offset),
                      ## this is deliberate, a better starting point
                      "cloglog" = glm.fit(x, y1, wt, family = binomial("probit"), offset = offset),
                      "cauchit" = glm.fit(x, y1, wt, family = binomial("cauchit"), offset = offset))
           if(!fit$converged)
               warning("attempt to find suitable starting values did not converge")
           coefs <- fit$coefficients
           if(any(is.na(coefs))) {
              warning("design appears to be rank-deficient, so dropping some coefs")
              keep <- names(coefs)[!is.na(coefs)]
              coefs <- coefs[keep]
 #          x <- x[, keep[-1], drop = FALSE]
              ## note: we keep the intercept
              nBeta <- ncol(x)
           }
           spacing <- logit((1:nInterval)/(nInterval+1)) # just a guess
           if(method != "logit") spacing <- spacing/1.7
           zetas <- -coefs[2] + spacing - spacing[q1]
           coefs[1] <- 0
           activePar <- c(FALSE, rep(TRUE, length(coef) - 1 + length(zetas) + 1))
                                         # intercept is fixed to 0
           sigma <- 1
        }
        else {
           ## not ordered: estimate OLS on interval middle points
           yMean <- numeric(length(yInt))
           if(nIntervalObs > 0) {
              means <- sapply(intervals, mean)
                            # we have to put a reasonable value to infinite intervals.
                            # Pick the average width of the interval and use it as the meanpoint
              widths <- sapply(intervals, function(x) x[2] - x[1])
              meanWidth <- mean(widths[!is.infinite(widths)])
              negInf <- is.infinite(means) & means < 0
              if(any(negInf)) {
                            # if none is true, sapply returns 'list()' and transforms means to a list
                 means[negInf] <- sapply(intervals[negInf], function(x) x[2] - meanWidth)
              }
              posInf <- is.infinite(means) & means > 0
              if(any(posInf)) {
                 means[posInf] <- sapply(intervals[posInf], function(x) x[1] + meanWidth)
              }
              yMean <- means[yInt]
           }
           yMean[!iIntervalObs] <- boundaries[boundaryInterval[!iIntervalObs]]
           yMean[is.infinite(yMean)] <- NA
                           # Inf will break lm ...
           fit <- lm(yMean ~ x - 1, na.action=na.action)
           xCoefs <- coef(fit)
           if(any(is.na(xCoefs))) {
              cat("Suggested initial values:\n")
              print(xCoefs)
              stop("NA in the initial values")
           }
           names(xCoefs) <- gsub("^x", "", names(xCoefs))
           sigma <- sqrt(var(residuals(fit)))
        }
        start[iBeta] <- xCoefs
        names(start)[iBeta] <- names(xCoefs)
        start[iBoundaries] <- boundaries
        names(start)[iBoundaries] <- names(boundaries)
        start[iSigma] <- sigma
        names(start)[iSigma] <- "sigma"
     }
     else {
        if(length(start) != iSigma)
           stop("'start' is of wrong length:\n",
                "The current model includes ", nBeta,
                " explanatory variables plus\n",
                length(iBeta), " interval boundaries ",
                "plus 1 disturbance standard deviation\n",
                "(", iSigma, " in total).\n",
                "However, 'start' is of length ",
                length(start))
     }
     if(print.level > 0) {
        cat("Initial values:\n")
        print(start)
     }
     if(!ordered) {
        ## Not ordered model: fix the fixed parameters
        activePar <- logical(length(start))
        activePar[iBeta] <- TRUE
        activePar[iBoundaries] <- FALSE
        activePar[iSigma] <- TRUE
     }
    ## compareDerivatives(loglik, gradlik, t0=start)
    ## stop()
    estimation <- maxLik(loglik, gradlik,
                         start=start,
                         method="BHHH", activePar=activePar,
                         control=list(iterlim=iterlim,
                         printLevel=print.level),
                         ...)
    res <- c(estimation,
             param=list(list(ordered=ordered,
                      boundaries=boundaries,
                      index=list(beta=iBeta, boundary=iBoundaries, std=iSigma),
                      intervalObs=iIntervalObsInf,
                           # include the Inf obs as these will be passed to
                           # the model.frame
                      df=nObs - sum(activePar),
                      nObs=nObs
                      )),
             call=cl,
             terms=mt,
             method=method,
             na.action=list(attr(mf, "na.action"))
             )
    class(res) <- c("intReg", class(estimation))
    return(res)
}
