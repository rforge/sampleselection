probit <- function(formula, ...) {
   ## Probit binary choice model.  Essentially a wrapper for "binaryCoice"
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
   ##
   ## we set the necessary function for binary choice gradient, and call it
   cdfLower <- pnorm
   cdfUpper <- function(x) pnorm(x, lower.tail=FALSE)
   logCdfLower <- function(x) pnorm(x, log.p=TRUE)
   logCdfUpper <- function(x) pnorm(x, lower.tail=FALSE, log.p=TRUE)
   pdf <- dnorm
   logPdf <- function(x) dnorm(x, log=TRUE)
   gradPdf <- function(x) -dnorm(x)*x
                           # these are the probit-specific distribution functions
   result <- binaryChoice(formula, ...,
                          cdfLower=cdfLower, cdfUpper=cdfUpper,
                          logCdfLower=logCdfLower, logCdfUpper=logCdfUpper,
                          pdf=pdf, logPdf=logPdf,
                          gradPdf=gradPdf)
   cl <- class(result)
   result <- c(result,
               family=list(binomial(link="probit"))
                           # NA action and the removed cases
               )
   class(result) <- c("probit", cl)
   result
}
