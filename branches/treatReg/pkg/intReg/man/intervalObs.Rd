\name{intervalObs}
\alias{intervalObs}
\alias{intervalObs.intReg}
\alias{minIntervalWidth}
\alias{minIntervalWidth,intReg-method}
\title{
  which observations are treated as interval/point observations
}
\description{
  \code{intervalObs} returns which observations are treated as interval or point
  observations
  when estimating
  \dQuote{intReg} objects.  \code{minIntervalWidth} returns the minimal
  width under which the intervals are treated as point observations.
}
\usage{
\method{intervalObs}{intReg}(object)
\S4method{minIntervalWidth}{intReg}(object, \dots)
}
\arguments{
  \item{object}{
    \dQuote{intReg} object
  }
  \item{\dots}{
    currently ignored
  }
}
\value{
  \code{intervalObs} returns a
  logical vector, which observations are treated as interval
  observations.  Others are treated as point observations.  Infinite
  intervals, such as \eqn{[a, \infty)}{[a, \code{Inf})}
  are treated as
  interval observations.

  \code{minIntervalWidth} returns the minimal interval width under which
  the observations are treated as point observations.  It is the same
  number as supplied by \code{minIntervalWidth} argument to
  \code{\link{intReg}}. 
}
\author{
  Ott Toomet
}
\keyword{methods}
\keyword{utilities}
