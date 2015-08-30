\name{intervalObs}
\alias{intervalObs}
\alias{intervalObs.intReg}
\title{
  which observations are treated as interval/point observations
}
\description{
  return which observations are treated as interval or point
  observations
  when estimating
  \dQuote{intReg} objects.
}
\usage{
\method{intervalObs}{intReg}(object)
}
\arguments{
  \item{object}{
    \dQuote{intReg} object
}
}
\value{
  Logical vector, which observations are treated as interval
  observations.  Others are treated as point observations.  Infinite
  intervals, such as \eqn{[a, \infty)}{[a, \code{Inf})} are treated as
  point observations.
}
\author{
  Ott Toomet
}
\keyword{methods}
\keyword{utilities}
