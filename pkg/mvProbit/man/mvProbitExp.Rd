\name{mvProbitExp}
\alias{mvProbitExp}
\title{Expectations from Multivariate Probit Models}
\description{
   Calculating expected outcomes from multivariate probit models.
}

\usage{
mvProbitExp( formula, coef, sigma, data,
   cond = FALSE, random.seed = 123, ... )
}

\arguments{
   \item{formula}{a one-sided or two-sided formula \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have the same regressors).}
   \item{coef}{a numeric vector of all coefficients.}
   \item{sigma}{the covariance/correlation matrix of the residuals
      (must be symmetric and have ones on its diagonal).}
   \item{data}{a \code{data.frame} containing the data.}
   \item{cond}{logical values indicating 
      whether conditional expectations (if \code{TRUE})
      or unconditional expectations (if \code{FALSE}, default)
      should be returned.}
   \item{random.seed}{an integer used to seed R's random number generator;
      this is to ensure replicability 
      when computing the (cumulative) multivariate normal distribution
      which is required to calculate the conditional expectations;
      defaults to 123.}
   \item{\dots}{additional arguments are passed
      to \code{\link[mvtnorm]{pmvnorm}}
      when calculating conditional expectations.}
}

\value{
   A data frame containing the expectations of the dependent variables.
}

\details{
   When calculating unconditional expectations,
   the left-hand side of argument \code{formula} is ignored.
   When calculating unconditional expectations
   and argument \code{formula} is a one-sided formula
   (i.e. only the right-hand side is specified),
   the conditional expectations are calculated based on the assumption
   that all other dependent variables are one.

   The \sQuote{state} (or \sQuote{seed}) of R's random number generator 
   is saved at the beginning of the \code{mvProbitExp} function 
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

\seealso{\code{\link[sampleSelection]{probit}},
   \code{\link[stats]{glm}}}

\examples{
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10

# generate explanatory variables
xData <- data.frame( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = as.numeric( rnorm( nObs ) > 0 ),
   x3 = rnorm( nObs ),
   x4 = rnorm( nObs ) )

# coefficients
beta <- c(  0.8,  1.2, -1.0,  1.4, -0.8,
           -0.6,  1.0,  0.6, -1.2, -1.6,
            0.5, -0.6, -0.7,  1.1,  1.2 )

# covariance matrix of error terms
library( "miscTools" )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# unconditional expectations of dependent variables
yExp <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = xData )
print( yExp )

# conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
yExpCond <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = beta,
   sigma = sigma, data = xData, cond = TRUE )
print( yExpCond )

# conditional expectations of dependent variables
# (assuming that all other dependent variables are zero)
xData$y1 <- 0
xData$y2 <- 0
xData$y3 <- 0
yExpCond0 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = beta, sigma = sigma, data = xData, cond = TRUE )
print( yExpCond0 )
}

\keyword{models}
\keyword{regression}
