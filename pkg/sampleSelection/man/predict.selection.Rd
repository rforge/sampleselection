\name{predict.selection}
\alias{predict.selection}

\title{Predict method for fitted sample selection models}
\description{
   Calculate predicted values for sample selection models
   fitted with function \code{\link{selection}}.
}

\usage{
   \method{predict}{selection}( object, newdata = NULL,
   part = ifelse( type \%in\% c( "unconditional", "conditional" ),
      "outcome", "selection" ),
   type = "unconditional", ... )
}

\arguments{
   \item{object}{a fitted object of class \code{selection}.}
   \item{newdata}{optionally, a data frame in which to look for variables
      with which to predict.
      If omitted, the fitted linear predictors or the fitted response values
      are returned.}
   \item{part}{character string indicating for which equation
      the predicted variables should be calculated:
      \code{"selection"} for the predicted values of the selection equation and
      \code{"outcome"} for the predicted values of the outcome equation;
      this argument is automatically chosen depending on the value
      of argument \code{type}.}
   \item{type}{if argument \code{type} is \code{"link"}
      and argument \code{part} is \code{"selection"},
      the linear predictors of the selection equation are returned;
      if argument \code{type} is \code{"response"}
      and argument \code{part} is \code{"selection"},
      the predicted probabilities of the selection equation are returned;
      if argument \code{type} is \code{"unconditional"}
      and argument \code{part} is \code{"outcome"},
      the unconditional expectations are returned,
      i.e. \code{E[y|X] = X \%*\% beta};
      if argument \code{type} is \code{"conditional"}
      and argument \code{part} is \code{"outcome"},
      the conditional expectations are returned,
      i.e. \code{E[y|X,Z,w=1] = X \%*\% beta
      + rho * sigma * dnorm( Z \%*\% gamma ) / pnorm( Z \%*\% gamma )}.}
   \item{\dots}{further arguments (currently ignored).}
}

\value{
   A numeric vector of the predicted values.
}

\author{
   Arne Henningsen and \sQuote{fg nu}
   (the code is partly based on the code posted by \sQuote{fg nu} at
   \url{http://stackoverflow.com/questions/14005788/predict-function-for-heckman-model})
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
predict( m, type = "conditional" ) 
predict( m, type = "link" ) 
predict( m, type = "response" ) 
predict( m, newdata = Mroz87[ 3:9, ] ) 
}

\keyword{methods}
