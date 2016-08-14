\name{intReg-package}
\alias{intReg-package}
\docType{package}
\title{
  Interval Regression 
}
\description{
  This package estimates interval regression.  Different disturbance
  distributions are supported and methods for predicting
  different outcomes are included.
}
\details{
  Interval regression is a linear regression model where the outcome
  variable is only observable in intervals.  Assume \eqn{y^*}{y*} is the
  latent outcome variable (link).  Instead of \eqn{y^*}{y*} we observe
  \deqn{y = (l, u) \quad\mbox{where}\quad l \le y^* < u.}{y = (l, y)
    where l <= y < y.}  \eqn{l} and \eqn{u} are the lower and upper
  boundaries, and we know that between lies between those.  Note that in
  case \eqn{l=y}, we have a point observation.

  If \eqn{y^*}{y*} is related to the explanatory variables \eqn{x} as in
  the ordinary linear model
  \deqn{y^* = \beta'x + \epsilon}{y* = \beta'x + \epsilon,} the
  parameters \eqn{\beta} can be estimated by Maximum Likelihood.

  Unlike the ordinary linear regression, there are more values one can
  predict in interval regression.  intReg supports predicting:
  \itemize{
    \item Link: \eqn{y^*}{y*}
    \item Conditional link: expected \eqn{y^*}{y*} value, given it falls
    in the interval \eqn{(l,u)}: \eqn{\mathrm{E} [y^*| l \le y^* < u]}{E[y*|l
      <= y* < u]}
  }
  
  Interval Regression shares similarities with several other censored
  linear model.  Probit (and logit) are essentially interval regression
  models with only two observed intervals, \eqn{(-\infty, 0)}{(-Inf, 0)} and
  \eqn{[0, \infty)}{[0, Inf)}.  Ordered regression is a form of interval
  regression where the interval boundaries are unobserved and must be
  estimated.  Finally, ordinary linear regression is an interval
  regression where all observations are point observations.

  The central function of the package is \code{\link{intReg}} that
  estimates the model.  It includes methods for \code{\link{coef}},
  \code{\link{summary}}, \code{\link{predict}}, and others.
}
\author{
Ott Toomet <otoomet@gmail.com>, with original source adapted from R core team.

Maintainer: Ott Toomet <otoomet@gmail.com>
}
%% \references{
%% }
\keyword{Statistics|regression}
\keyword{Statistics|multivariate}
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{
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
}
