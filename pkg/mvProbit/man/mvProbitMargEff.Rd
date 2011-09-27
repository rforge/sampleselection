\name{mvProbitMargEff}
\alias{mvProbitMargEff}
\title{Marginal Effects based on Multivariate Probit Models}
\description{
   Calculating marginal effects based on multivariate probit models.
}

\usage{
mvProbitMargEff( formula, coef, sigma, data,
   cond = FALSE, algorithm = GenzBretz(), eps = 1e-06, random.seed = 123,
   ... )
}

\arguments{
   \item{formula}{a one-sided or two-sided \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have the same regressors).}
   \item{coef}{a numeric vector of all coefficients.}
   \item{sigma}{the covariance/correlation matrix of the residuals
      (must be symmetric and have ones on its diagonal).}
   \item{data}{a \code{data.frame} containing the data.}
   \item{cond}{logical value indicating 
      whether marginal effects on conditional expectations (if \code{TRUE})
      or marginal effects on unconditional expectations (if \code{FALSE}, default)
      should be returned.}
   \item{algorithm}{algorithm for computing integrals 
      of the multivariate normal distribution,
      either \code{GenzBretz()}, \code{Miwa()}, or \code{TVPACK()}
      (see documentation of \code{\link[mvtnorm]{pmvnorm}}).}
   \item{eps}{numeric, the difference between the two values 
      of each (numerical) explanatory
      variable that is used for the numerical differentiation.}
   \item{random.seed}{an integer used to seed R's random number generator;
      this is to ensure replicability 
      when computing (cumulative) probabilities of the multivariate normal distribution
      which is required to calculate the conditional expectations;
      \code{set.seed( random.seed )} is called each time before
      a (cumulative) probability of the multivariate normal distribution 
      is computed;
      defaults to 123.}
   \item{\dots}{additional arguments are passed
      to \code{\link[mvtnorm]{pmvnorm}}
      when calculating conditional expectations.}
}

\value{
   A data frame containing the marginal effects of the explanatory variables
   on the expectations of the dependent variables.
}

\details{
   When calculating marginal effects on unconditional expectations,
   the left-hand side of argument \code{formula} is ignored.
   When calculating marginal effects on unconditional expectations
   and argument \code{formula} is a one-sided formula
   (i.e. only the right-hand side is specified),
   the marginal effects on the conditional expectations 
   are calculated based on the assumption
   that all other dependent variables are one.

   If a dependent variable contains only zeros and ones 
   or only \code{TRUE} and \code{FALSE},
   argument \code{eps} is ignored (for this variable) 
   and the effect between this variable being zero and one
   is returned.

   The \sQuote{state} (or \sQuote{seed}) of R's random number generator 
   is saved at the beginning of the \code{mvProbitMargEff} function 
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

\seealso{\code{\link{mvProbit}},
   \code{\link{mvProbitExp}},
   \code{\link{mvProbitLogLik}},
   \code{\link[sampleSelection]{probit}},
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

# marginal effects on unconditional expectations of dependent variables
margEffUnc <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = xData )
print( margEffUnc )

# marginal effects on conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
margEffCond <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = beta,
   sigma = sigma, data = xData, cond = TRUE )
print( margEffCond )

# marginal effects on conditional expectations of dependent variables
# (assuming that all other dependent variables are zero)
xData$y1 <- 0
xData$y2 <- 0
xData$y3 <- 0
margEffCond0 <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = beta, sigma = sigma, data = xData, cond = TRUE )
print( margEffCond0 )
}

\keyword{models}
\keyword{regression}
