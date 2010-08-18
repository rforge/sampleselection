\name{tobit}
\alias{tobit}
\title{Censored Regression (Tobit) Model}
\description{
   Fitting a model with a censored dependent variable.
}

\usage{
tobit( formula, left = 0, right = Inf, data = sys.frame( sys.parent()),
   start = NULL, nGHQ = 4, ... )
}

\arguments{
   \item{formula}{an object of class \code{"formula"}:
      a symbolic description of the model to be fitted.}
   \item{left}{left limit for the censored dependent variable;
      if set to \code{-Inf}, the dependent variable is assumed to be
      not left-censored; defaults to zero (classical Tobit model).}
   \item{right}{right limit for the censored dependent variable;
      if set to \code{Inf}, the dependent variable is assumed to be
      not right-censored; defaults to \code{Inf} (classical Tobit model).}
   \item{data}{an optional data frame.
      If argument \code{data} is of class \code{"pdata.frame"},
      a panel-tobit model is estimated.}
   \item{start}{an optional vector of initial parameters for the ML estimation
      (intercept, slope parameters, and logarithm of the standard deviation
      of the disturbance term);
      if \code{start} is not specified, initial values are taken from an OLS
      estimation.}
   \item{nGHQ}{number of points used in the Gauss-Hermite quadrature,
      which is used to compute the log-likelihood value in case of the
      random effects model for panel data.} 
  \item{\dots}{additional arguments are passed to
      \code{\link[maxLik]{maxLik}}.}
}

\value{
   \code{tobit} returns an object of class \code{"tobit"}
   inheriting from class \code{"maxLik"}.
}

\details{
   The model is estimated by Maximum Likelihood (ML)
   assuming a Gaussian (normal) distribution of the error term.
   The maximization of the likelihood function is done
   by function \code{\link[maxLik]{maxLik}} of the \pkg{maxLik} package.
   An additional argument \code{method} can be used to specify
   the optimization method used by \code{\link[maxLik]{maxLik}},
   e.g.\ \code{"Newton-Raphson"}, \code{"BHHH"}, \code{"BFGS"},
   \code{"SANN"} (for simulated annealing), or
   \code{"NM"} (for Nelder-Mead).
}

\references{
   Greene, W.H. (2008):
   \emph{Econometric Analysis}, Sixth Edition, Prentice Hall, p. 871-875.

   Kleiber, C. and Zeileis, A. (2008):
   \emph{Applied Econometrics with R}, Springer, p. 141-143.

   Tobin, J. (1958):
   Estimation of Relationships for Limited Dependent Variables.
   \emph{Econometrica} 26, p. 24-36.
}

\note{
  The 2-step estimate of 'rho' may be outside of the
  \eqn{[-1,1]}{[-1,1]} interval.  In that case the standard errors of
  invMillsRatio may be meaningless.
  }

\author{Arne Henningsen}

\seealso{\code{\link[AER]{tobit}}, \code{\link[sampleSelection]{selection}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
tobitResult <- tobit( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
summary( tobitResult )

## Kleiber & Zeileis ( 2008 ), page 143
tobitResultBoth <- tobit( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, right = 4 )
summary( tobitResultBoth )
}

\keyword{models}
\keyword{regression}
