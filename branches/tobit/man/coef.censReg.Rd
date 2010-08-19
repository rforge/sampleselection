\name{coef.censReg}
\alias{coef.censReg}

\title{Extract Coefficients from Censored Regression Models}

\description{
   This function extracts coefficients
   from censored regression models
}

\usage{
\method{coef}{censReg}( object, logSigma = TRUE, ... )
}

\arguments{
   \item{object}{object of class \code{censReg}
      (returned by \code{\link{censReg}}).}
   \item{logSigma}{logical value indicating whether the variance(s)
      should be returned logarithmized.}
   \item{\dots}{currently not used.}
}
\value{
   \code{coef.censReg} returns a vector of the estimated coefficients.
}

\note{
   When the censored regression model is estimated,
   the log-likelihood function is maximized with respect
   to the coefficients and the \emph{logarithm(s)}
   of the variance(s).
   Hence, if argument \code{logSigma} is \code{FALSE},
   the variance(s) is/are calculated by applying the exponential function
   to the logarithmized variance(s).
}

\author{Arne Henningsen}

\seealso{\code{\link{censReg}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
coef( estResult )
coef( estResult, logSigma = TRUE )
}

\keyword{methods}
