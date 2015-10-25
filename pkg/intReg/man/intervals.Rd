\name{intervals}
\alias{intervals}
\alias{intervals.intReg}
\alias{intervals.list}
\alias{interva.data.frame}
\title{
  return interval boundaries for each observation
}
\description{
  find interval boundaries for estimated interval regression models.
}
\usage{
intervals(object)
}
\arguments{
  \item{object}{
    'intReg' object, or a list or data frame.  List (and data frame)
    must contain either: a component \sQuote{yInt}, factor in the form
    \code{(LB, UB)} (brackets instead of parentheses are allowed) where
    \sQuote{LB} and \sQuote{UB} are numeric values for lower and upper bound; or
    numeric components \sQuote{LB} and \sQuote{UB}.  Output from
    \code{\link{cut}} is suitable for the first form.  If both forms are
    present, the former is used.
  }
}
\value{
  a matrix with columns \sQuote{LB} and \sQuote{UB} for lower and upper
  bound respectively.
}
\author{
  Ott Toomet
}
\keyword{methods}
\keyword{utilities}
