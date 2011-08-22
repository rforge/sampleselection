\name{mvProbitLogLik}
\alias{mvProbitLogLik}
\title{Log Likelihood Values for Multivariate Probit Models}
\description{
   Calculating log likelihood values for multivariate probit models.
}

\usage{
mvProbitLogLik( formula, coef, sigma, data, random.seed = 123, ... )
}

\arguments{
   \item{formula}{a \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have the same regressors).}
   \item{coef}{a numeric vector of all coefficients.}
   \item{sigma}{the covariance/correlation matrix of the residuals
      (must be symmetric and have ones on its diagonal).}
   \item{data}{a \code{data.frame} containing the data.}
   \item{random.seed}{an integer used to seed R's random number generator;
      this is to ensure replicability 
      when computing the (cumulative) multivariate normal distribution
      which is required to calculate the log likelihood values;
      defaults to 123.}
   \item{\dots}{additional arguments are passed
      to \code{\link[mvtnorm]{pmvnorm}}
      when calculating conditional expectations.}
}

\value{
   A vector containing the log likelihood values for each observation.
}

\details{
   The \sQuote{state} (or \sQuote{seed}) of R's random number generator 
   is saved at the beginning of the \code{mvProbitLogLik} function 
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

\seealso{\code{\link{mvProbitExp}},
   \code{\link{mvProbitMargEff}},
   \code{\link[sampleSelection]{probit}},
   \code{\link[stats]{glm}}}

\examples{
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = as.numeric( rnorm( nObs ) > 0 ),
   x3 = rnorm( nObs ),
   x4 = rnorm( nObs ) )

# coefficients
beta <- cbind( c(  0.8,  1.2, -1.0,  1.4, -0.8 ),
               c( -0.6,  1.0,  0.6, -1.2, -1.6 ),
               c(  0.5, -0.6, -0.7,  1.1,  1.2 ) )

# covariance matrix of error terms
library( "miscTools" )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat \%*\% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )

# log likelihood values
myData <- as.data.frame( cbind( xMat, yMat ) )
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = myData )
print( logLikVal )
}

\keyword{models}
\keyword{regression}
