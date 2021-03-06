\name{nObs.intReg}
\alias{nObs.intReg}
\title{
  Number of observations.
}
\description{
  Extracts the number of observations in the intReg model
}
\usage{
\method{nObs}{intReg}(x, ...)
}
\arguments{
  \item{x}{
    object of class \sQuote{intReg}, estimated interval regression model
  }
  \item{\dots}{other parameters for different methods}
}
\value{
  Integer, number of observations
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
Kakadu <- Kakadu[sample(nrow(Kakadu), 500),]
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
print(nObs(m))
}
\keyword{methods}
