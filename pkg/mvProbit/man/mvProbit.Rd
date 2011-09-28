\name{mvProbit}
\alias{mvProbit}
\title{Estimation of Multivariate Probit Models}
\description{
   Estimating multivariate probit models by the maximum likelihood method.
   
   WARNING: this function is experimental and extremely
   (perhaps even unusably) slow!
}

\usage{
mvProbit( formula, coef = NULL, sigma = NULL, data, method = "BHHH", finalHessian = "BHHH", 
   algorithm = GenzBretz(), nGHK = 1000, oneSidedGrad = FALSE, eps = 1e-6, 
   random.seed = 123, ... )
}

\arguments{
   \item{formula}{a \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have the same regressors).}
   \item{coef}{an optional numeric vector specifying the starting values 
      of all coefficients; 
      if this argument is not specified,
      coefficients estimated by univariate probit models 
      are used as starting values for the coefficients.}
   \item{sigma}{optional starting values for the covariance/correlation matrix 
      of the residuals (must be symmetric and have ones on its diagonal);
      if this argument is not specified,
      the correlation matrix of the \sQuote{response} residuals, 
      i.e. y - pnorm( X' beta ),
      is used as starting values for sigma.}
   \item{data}{a \code{data.frame} containing the data.}
   \item{method}{maximisation method / algorithm 
      (see \code{\link[maxLik]{maxLik}}).}
   \item{finalHessian}{Calculation of the final Hessian:
      either \code{FALSE} (no calculation of Hessian), 
      \code{TRUE} (finite-distance calculation of Hessian), or 
      \code{"BHHH"} (calculation based on information equality approach
         and finite-distance gradients, the default).}
   \item{algorithm}{algorithm for computing integrals 
      of the multivariate normal distribution,
      either function \code{GenzBretz()}, \code{Miwa()}, or \code{TVPACK()}
      (see documentation of \code{\link[mvtnorm]{pmvnorm}})
      or character string \code{"GHK"} 
      (see documentation of \code{\link[bayesm]{ghkvec}}).}
   \item{nGHK}{numeric value specifying the number of simulation draws 
      of the GHK algorithm for computing integrals 
      of the multivariate normal distribution.}
   \item{oneSidedGrad}{logical. If \code{TRUE}, 
      the gradients of the log-likelihood function
      with respect to the estimated parameters 
      are obtained by one-sided numeric finit-difference differentiation,
      which is faster but less precise 
      than two-sided numeric finit-difference differentiation.
      This argument should be set to \code{FALSE}
      if an optimisation algorithm is used that is not based on gradients.}
   \item{eps}{numeric. The step size for the one-sided numeric
      finit-distance differentiation.
      Unfortunately, it is currently not possible to set the step size
      for the two-sided numeric finit-distance differentiation.}
   \item{random.seed}{an integer used to seed R's random number generator;
      this is to ensure replicability 
      when computing (cumulative) probabilities of the multivariate normal distribution
      which is required to calculate the log likelihood values;
      \code{set.seed( random.seed )} is called each time before
      a (cumulative) probability of the multivariate normal distribution 
      is computed;
      defaults to 123.}
   \item{\dots}{additional arguments are passed
      to \code{\link[maxLik]{maxLik}} and \code{\link[mvtnorm]{pmvnorm}}.}
}

\value{
   The object returned by \code{\link[maxLik]{maxLik}}.
}

\details{
   The \sQuote{state} (or \sQuote{seed}) of R's random number generator 
   is saved at the beginning of the \code{mvProbit} function 
   and restored at the end of this function 
   so that this function does \emph{not} affect the generation 
   of random numbers outside this function
   although the random seed is set to argument \code{random.seed}
   and the calculation of the (cumulative) multivariate normal distribution
   uses random numbers.
}

\references{
   Greene, W.H. (1996): 
   \emph{Marginal Effects in the Bivariate Probit Model},
   NYU Working Paper No. EC-96-11. 
   Available at \url{http://ssrn.com/abstract=1293106}.
}

\author{Arne Henningsen}

\seealso{\code{\link{mvProbitLogLik}},
   \code{\link{mvProbitExp}},
   \code{\link{mvProbitMargEff}},
   \code{\link[sampleSelection]{probit}},
   \code{\link[stats]{glm}}}

\examples{
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 25

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = rnorm( nObs ) )

# coefficients
beta <- cbind( c(  0.8,  1.2, -0.8 ),
               c( -0.6,  1.0, -1.6 ),
               c(  0.5, -0.6,  1.2 ) )

# covariance matrix of error terms
library( "miscTools" )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat \%*\% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )

# estimation with the BFGS algorithm
estResultBFGS <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5 )
summary( estResultBFGS )
}

\keyword{models}
\keyword{regression}
