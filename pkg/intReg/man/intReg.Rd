\name{intReg}
\alias{intReg}
\alias{print.intReg}
\alias{print.summary.intReg}
\alias{summary.intReg}
\title{Interval Regression}
\description{
  This function estimates interval regression using either common
  intervals for all observations, or observation-specific intervals.
  Normal, logistic, log-log, and Cauchy disturbances are supported.
}
\usage{
intReg(formula, start, boundaries, ...,
contrasts = NULL, Hess = FALSE, model = TRUE,
method = c("probit", "logistic", "cloglog", "cauchit", "model.frame"), print.level = 0,
data, subset, weights, na.action,
iterlim=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class '"formula"' (or one that can be coerced to
          that class): a symbolic description of the model to be
          fitted.  The left-hand-side variable must be either a factor,
	  describing the intervals where the observations fall to, or a matrix of
	  two columns, describing the interval boundaries for each
	  observation.  See details below, see
	  also \code{\link{lm}} for explanation how to write formulas.}
  \item{data}{an optional data frame, list or environment (or object
          coercible by 'as.data.frame' to a data frame) containing the
          variables in the model.  If not found in 'data', the
          variables are taken from 'environment(formula)', typically
          the environment from which \sQuote{intReg} is called.}
  \item{weights}{an optional vector of weights to be used in the fitting
          process.  Should be 'NULL' or a numeric vector. If non-NULL,
          weighted likelihood is maximized (that
          is, maximizing \sQuote{sum(w*loglik)}).}
  \item{start}{Initial values for the optimization algorithm.  if
    \sQuote{NULL}, initial values are calculated based on interval
    means.  Note that \code{intReg} expects the full-length inital
    values, including the parameter boundaries and the standard
    deviation of the error term.  See the example below.
  }
  \item{boundaries}{Boundaries for intervals.  See details.}
  \item{\dots}{further arguments to \code{\link[maxLik]{maxLik}}}
  \item{subset}{an optional logical vector specifying a subset of observations to be
          used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data
          contain \sQuote{NA}s.  The default is set by the \sQuote{na.action} setting
          of sQuote{options}, and is \sQuote{na.fail} if that is unset.  The
          \dQuote{factory-fresh} default is \sQuote{na.omit}.  Another possible value
          is \sQuote{NULL}, no action.  Value \sQuote{na.exclude} can be useful.}
  \item{contrasts}{an optional list. See the \sQuote{contrasts.arg} of
          \code{\link{model.matrix.default}}.}
  \item{Hess}{Should Hessian of the model be returned}
  \item{model}{logical for whether the model matrix should be returned.}
  \item{method}{character, distribution of disturbances or \sQuote{model.frame}.
    \sQuote{probit}, \sQuote{logistic}, \sQuote{cloglog} and
    \sQuote{cauchit} assume these disturbance distributions.
    \sQuote{model.frame} returns the model frame.  The default value
    \sQuote{probit} assumes normal distribution.}
  \item{print.level}{output of run-time information: higher level prints
    more.}
  \item{iterlim}{maximum number of optimization iterations}
}
\details{
  Interval regression is a form of linear regression where only
  intervals (i.e numeric upper and lower bounds) where the observations fall are
  visible of the otherwise continuous outcome
  variable.  The current implementation assumes known distribution of
  the error term (default is normal) and estimates the model using
  Maximum Likelihood.

  The intervals may be specified in two ways: either common intervals
  for all the observations by using argument 'boundaries', or by
  specifying the response as Nx2 matrix, where columns correspond to the
  lower- and upper bound for the individual observations.  If you
  specify the common boundaries, you may want to give names to the
  components of 'boundaries'.  Otherwise they will be names
  \sQuote{Boundary 1}, \sQuote{Boundary 2} etc.  These names will appear to the estimation
  results.  For observations-specific boundaries, the names are
  generated automatically.
}
\value{
  Object of class \sQuote{intReg} which inherits from class \sQuote{maxLik}.

  There are several methods, including \code{\link{summary}} and
  \code{\link{coef}}, partly inherited from "maxLik" class.
  
  Note that the boundaries are passed as fixed parameters to the
  \code{maxLik} estimation routine and hence returned as fixed estimates
  with standard errors set to zero.
}
\author{Ott Toomet \email{otoomet@gmail.com}}
\note{A lot of code is borrowed from \code{\link[MASS]{polr}} from the package MASS
}
\seealso{\code{\link[MASS]{polr}}}
\examples{
## Example of observation-specific boundaries
## Estimate the willingness to pay for the Kakadu National Park
## Data given in intervals -- 'lower' for lower bound and 'upper' for upper bound.
## Note that dichotomous-coice answers are already coded to 'lower' and 'upper'
data(Kakadu, package="Ecdat")
set.seed(1)
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
print(summary(m))
##
## Example of setting initial values
##
st <- coef(m, boundaries=TRUE)
st[1:19] <- 1     # set all coefficients to 1
st["sigma"] <- 1  # set standard deviation to 1
m <- intReg(y ~ sex + log(income) + age + schooling + 
              recparks + jobs + lowrisk + wildlife + future + aboriginal + finben +
              mineparks + moreparks + gov +
              envcon + vparks + tvenv + major,
              start=st,
              data=Kakadu)
## Note: the results will be the same as above
##
## Example of common intervals for all the observations
##
library(Ecdat)
data(Bwages)
## calculate an ordinary Mincer-style wage regression.  First by OLS and
## thereafter cut the wage to intervals and estimate with 'intReg'
## Note: gross hourly wage rate in EUR
ols <- lm(log(wage) ~ factor(educ) + poly(exper, 2), data=Bwages)
cat("OLS estimate:\n")
print(summary(ols))
## Now we censor the wages to intervals
intervals <- c(0, 5, 10, 15, 25, Inf)
salary <- cut(Bwages$wage, intervals)
int <- intReg(salary ~ factor(educ) + poly(exper, 2), data=Bwages, boundaries=log(intervals))
## Note: use logs for the intervals in Euros.  We do not have to
## transform salaris to log form as this does not change the intervals.
## Ignore any warnings
cat("Interval regression:\n")
print(summary(int))
}
\keyword{models}
\keyword{regression}
