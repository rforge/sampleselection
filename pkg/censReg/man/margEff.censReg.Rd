\name{margEff.censReg}
\alias{margEff.censReg}
\alias{summary.margEff.censReg}

\title{Marginal Effects in Censored Regression Models}

\description{
   The \code{margEff} method computes the marginal effects of the explanatory variables
   on the expected value of the dependent variable evaluated at arbbitrary
   values of the explanatory variables, or at the mean values of the explanatory variables as the default.
   Please note that this functionality is currently not available
   for panel data models.
}

\usage{
\method{margEff}{censReg}( object, calcVCov = TRUE, returnJacobian = FALSE, values, ... )
# margEff(object)
# margEff(object, values="default2")
# margEff(object, values=c(...))

\method{summary}{margEff.censReg}( object, ... )
}

\arguments{
   \item{object}{argument \code{object} of the \code{margEff} method
      must be an object of class \code{"censReg"}
      (returned by \code{\link{censReg}});
      argument \code{object} of the \code{summary} method
      must be an object of class \code{"margEff.censReg"}
      (returned by \code{margEff.censReg}).}
   \item{calcVCov}{logical. If \code{TRUE}, 
     the approximate variance covariance matrices of the marginal effects 
     is calculated and returned as an attribute (see below).}
   \item{returnJacobian}{logical. If \code{TRUE}, 
      the Jacobian of the marginal effects with respect to the coefficients 
      is returned as an attribute (see below).}
   \item{values}{at what values of the independent variables marginal effects are calculated (see Value).} 
   \item{\dots}{currently not used.}
}

\value{
   \code{margEff.censReg} returns an object of class \code{"margEff.censReg"},
   which is a vector of the marginal effects of the explanatory variables. With  missing \code{values} it is calculated on the expected value of the dependent variable evaluated at the mean values of the explanatory variables. With 
   \code{values="default2"} means of continuos varibles and levels of factors coded with zeros are used. For factors given as functions of continuous variables such as \code{factor(x>0)} the effects are calculated at meanings of factors corresponding to means of continuous variables, e.g. \code{factor(mean(x)>0)}. Arbitrary values of the independent variables should be given as vectors. The values should have the same format as in the data for \code{censReg}. The order of variables should be as in \code{all.vars(object$call$formula)}.
   The returned object has an attribute \code{df.residual}, 
   which is equal to the degrees of freedom of the residuals. 
   
   If argument \code{calcVCov} is \code{TRUE}, 
   the object returned by \code{margEff.censReg} has an attribute \code{vcov}, 
   which is a the approximate variance covariance matrices 
   of the marginal effects calculated
   with the Delta method.
   
   If argument \code{returnJacobian} is \code{TRUE}, 
   the object returned by \code{margEff.censReg}has an attribute \code{jacobian}, 
   which is the Jacobian of the marginal effects 
   with respect to the coefficients.

   \code{summary.margEff.censReg} returns 
   an object of class \code{"summary.margEff.censReg"},
   which is a matrix with four columns,
   where the first column contains the marginal effects,
   the second column contains the standard errors of the marginal effects,
   the third column contains the corresponding t-values,
   and the fourth columns contains the corresponding P values.
}

\author{Arne Henningsen, Oleh Komashko}

\seealso{\code{\link{censReg}}, \code{\link{coef.censReg}}, 
   and \code{\link{summary.censReg}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
margEff( estResult )
margEff( estResult , values = "default2")
margEff( estResult , values = c(30,5,1,1,1) )
summary( margEff( estResult ) )
}

\keyword{methods}
