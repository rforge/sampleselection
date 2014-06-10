\name{predict.selection}
\alias{predict.selection}

\title{Predict method for fitted sample selection models}
\description{
   Calculate predicted values for sample selection models
   fitted with function \code{\link{selection}}.
}

\usage{
   \method{predict}{selection}( object, newdata = NULL, type = "uncond", ... )
}

\arguments{
   \item{object}{a fitted object of class \code{selection}.}
   \item{newdata}{optionally, a data frame in which to look for variables
      with which to predict.
      If omitted, the fitted linear predictors or the fitted response values
      are returned.}
   \item{type}{the type of prediction,
      \code{"uncond"} (the default), \code{"cond"},
      \code{"link"}, or \code{"prob"}.}
   \item{\dots}{further arguments (currently ignored).}
}

\value{
   A numeric vector of the predicted values.
}

\author{
   Arne Henningsen and \sQuote{fg nu}
   (the code is partly based on the code posted by \sQuote{fg nu} at
   \url{http://stackoverflow.com/questions/14005788/predict-function-for-heckman-model}
}

\seealso{
   \code{\link{selection}}, \code{\link[stats]{predict}},
   \code{\link{predict.probit}}, \code{\link{residuals.selection}},
   and \code{\link{selection-methods}}.
}

\examples{
## Greene( 2003 ): example 22.8, page 786
data( Mroz87 )
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )

# ML estimation
m <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
   
predict( m )
predict( m, type = "cond" ) 
predict( m, type = "link" ) 
predict( m, type = "prob" ) 
predict( m, newdata = Mroz87[ 3:9, ] ) 
}

\keyword{methods}
