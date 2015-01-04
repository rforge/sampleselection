\name{coef.intReg}
\alias{coef.intReg}
\alias{print.coef.intReg}
\alias{print.summary.intReg}
\alias{summary.intReg}
\alias{stdEr.intReg}
\alias{vcov.intReg}
\title{
  Extract (only informative) coefficients, standard errors, and
  variance-covariance matrix from \sQuote{intReg} model.
}
\description{
  Internally, the interval boundaries are treated as fixed parameters.
  These methods extract only the informative ones (by default),
  i.e. they do not return the interval boundaries.
}
\usage{
\method{coef}{intReg}(object, boundaries = FALSE, ...)
\method{stdEr}{intReg}(x, boundaries = FALSE, ...)
\method{summary}{intReg}(object, boundaries = FALSE, ...)
\method{vcov}{intReg}(object, boundaries = FALSE, ...)
}
\arguments{
  \item{object,x}{
    object of class \sQuote{intReg}, estimated interval regression model
}
  \item{boundaries}{
    logical, whether to return (fixed) interval boundary parameters
}
  \item{\dots}{
    arguments for other methods
}
}
\value{
  a named numeric vector or matrix, for the estimated coefficients and
  standard errors, and variance-covariance matrix respectively.
}
\author{
  Ott Toomet \email{otoomet@gmail.com}
}
\seealso{
  \code{\link{summary.intReg}} which provides related functionality.
}
\examples{
## Example of observation-specific boundaries
## Estimate the willingness to pay for the Kakadu National Park
## Data given in intervals -- 'lower' for lower bound and 'upper' for upper bound.
## Note that dichotomous-coice answers are already coded to 'lower' and 'upper'
data(Kakadu, package="Ecdat")
set.seed(1)
Kakadu <- Kakadu[sample(nrow(Kakadu), 400),]
                        # subsample to speed up the estimation
## Estimate in log form, change 999 to Inf
lb <- log(Kakadu$lower)
ub <- Kakadu$upper
ub[ub > 998] <- Inf
ub <- log(ub)
y <- cbind(lb, ub)
m <- intReg(y ~ sex + log(income) + age + schooling + 
              recparks + jobs + lowrisk + wildlife + future + aboriginal + finben +
              mineparks + moreparks + gov +
              envcon + vparks + tvenv + major, data=Kakadu)
## You may want to compare the results to Werner (1999),
## Journal of Business and Economics Statistics 17(4), pp 479-486
print(coef(m))
print(coef(m, boundaries=TRUE))
print(nObs(m))
}
\keyword{methods}
