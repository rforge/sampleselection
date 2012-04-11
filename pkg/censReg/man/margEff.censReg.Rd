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
\method{margEff}{censReg}( object, calcVCov = TRUE, returnJacobian = FALSE, ... )
}

\arguments{
   \item{object}{object of class \code{"censReg"}
      (returned by \code{\link{censReg}}).}
   \item{calcVCov}{logical. If \code{TRUE}, 
     the approximate variance covariance matrices of the marginal effects 
     is calculated and returned as an attribute (see below).}
   \item{returnJacobian}{logical. If \code{TRUE}, 
      the Jacobian of the marginal effects with respect to the coefficients 
      is returned as an attribute (see below).}
   \item{\dots}{currently not used.}
}

\value{
   \code{margEff.censReg} returns 
   the vector of the marginal effects of the explanatory variables
   on the expected value of the dependent variable calculated
   at the mean values of the explanatory variables.
   The returned object has an attribute \code{df.residual}, 
   which is equal to the degrees of freedom of the residuals. 
   
   If argument \code{calcVCov} is \code{TRUE}, 
   the returned object has an attribute \code{vcov}, 
   which is a the approximate variance covariance matrices 
   of the marginal effects calculated
   with the Delta method.
   
   If argument \code{returnJacobian} is \code{TRUE}, 
   the returned object has an attribute \code{jacobian}, 
   which is the Jacobian of the marginal effects 
   with respect to the coefficients.
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
