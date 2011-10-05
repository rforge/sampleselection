\name{mvProbitExp}
\alias{mvProbitExp}
\title{Expectations from Multivariate Probit Models}
\description{
   Calculating expected outcomes from multivariate probit models.
}

\usage{
mvProbitExp( formula, coef, sigma = NULL, data,
   cond = FALSE, algorithm = GenzBretz(), nGHK = 1000, random.seed = 123, ... )
}

\arguments{
   \item{formula}{a one-sided or two-sided \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have the same regressors).}
   \item{coef}{a numeric vector of all coefficients;
      if argument \code{sigma} is not specified,
      this vector must also include the correlation coefficients;
      the order of elements is explained in the section \dQuote{details}.}
   \item{sigma}{optional argument for specifying
      the covariance/correlation matrix of the residuals
      (must be symmetric and have ones on its diagonal);
      if this argument is not specified, 
      the correlation coefficients must be specified by argument \code{coef}.}
   \item{data}{a \code{data.frame} containing the data.}
   \item{cond}{logical value indicating 
      whether conditional expectations (if \code{TRUE})
      or unconditional expectations (if \code{FALSE}, default)
      should be returned.}
   \item{algorithm}{algorithm for computing integrals 
      of the multivariate normal distribution,
      either function \code{GenzBretz()}, \code{Miwa()}, or \code{TVPACK()}
      (see documentation of \code{\link[mvtnorm]{pmvnorm}})
      or character string \code{"GHK"} 
      (see documentation of \code{\link[bayesm]{ghkvec}}).}
   \item{nGHK}{numeric value specifying the number of simulation draws 
      of the GHK algorithm for computing integrals 
      of the multivariate normal distribution.}
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

   If the model has \eqn{n} dependent variables (equations)
   and \eqn{k} regressors in each equation,
   the order of the coefficients in argument \code{coef} must be as follows:
   \eqn{b_{1,1}}, \ldots, \eqn{b_{1,k}}, 
   \eqn{b_{2,1}}, \ldots, \eqn{b_{2,k}}, \ldots,
   \eqn{b_{n,1}}, \ldots, \eqn{b_{n,k}},
   where \eqn{b_{i,j}} is the coefficient
   of the \eqn{j}th regressor in the \eqn{i}th equation.
   If argument \code{sigma} is not specified,
   argument \code{coef} must additionally include following elements:
   \eqn{R_{1,2}}, \eqn{R_{1,3}}, \eqn{R_{1,4}}, \ldots, \eqn{R_{1,n}},
   \eqn{R_{2,3}}, \eqn{R_{2,4}}, \ldots, \eqn{R_{2,n}}, \ldots,
   \eqn{R_{n-1,n}},
   where \eqn{R_{i,j}} is the correlation coefficient corresponding to
   the \eqn{i}th and \eqn{j}th equation.

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

\seealso{\code{\link{mvProbit}},
   \code{\link{mvProbitMargEff}},
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
