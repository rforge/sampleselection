\name{summary.censReg}
\alias{summary.censReg}

\title{Summary Results of Censored Regression Models}

\description{
   This function prepares summary results
   for censored regression models.
}

\usage{
\method{summary}{censReg}( object, ... )
}

\arguments{
   \item{object}{object of class \code{"censReg"}
      (returned by \code{\link{censReg}}).}
   \item{\dots}{currently not used.}
}
\value{
   \code{summary.censReg} returns an object of class \code{"summary.censReg"}
   inheriting from class \code{"summary.maxLik"}.
   The returned object contains the same components as objects
   returned by \code{\link[maxLik]{summary.maxLik}} and additionally
   the following components:
   \item{call}{the matched call.}
   \item{nObs}{a vector containing 4 integer values:
      the total number of observations,
      the number of left-censored observations,
      the number of uncensored observations, and
      the number of right-censored observations.}
}

\author{Arne Henningsen}

\seealso{\code{\link{censReg}} and \code{\link{coef.censReg}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
summary( estResult )
}

\keyword{methods}
