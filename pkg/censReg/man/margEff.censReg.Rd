\name{margEff.censReg}
\alias{margEff.censReg}

\title{Marginal Effects in Censored Regression Models}

\description{
   This function computes the marginal effects of the explanatory variables
   on the expected value of the dependent variable calculated
   at the mean values of the explanatory variables.
   Please note that this functionality is currently not available
   for panel data models.
}

\usage{
\method{margEff}{censReg}( object, ... )
}

\arguments{
   \item{object}{object of class \code{"censReg"}
      (returned by \code{\link{censReg}}).}
   \item{\dots}{currently not used.}
}
\value{
   \code{margEff.censReg} returns 
   the vector of the marginal effects of the explanatory variables
   on the expected value of the dependent variable calculated
   at the mean values of the explanatory variables.
}

\author{Arne Henningsen}

\seealso{\code{\link{censReg}}, \code{\link{coef.censReg}}, 
   and \code{\link{summary.censReg}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
margEff( estResult )
}

\keyword{methods}
