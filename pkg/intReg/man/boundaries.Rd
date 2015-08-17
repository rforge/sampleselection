\name{boundaries}
\alias{boundaries}
\title{
  Return interval boundaries for intReg models 
}
\description{
  Return a vector of all interval boundaries in increasing order.  This
  also includes eventual \code{-Inf} and \code{Inf} if used.  In case of
  observation-specific boundaries, no
  distinction is made if a boundary is used as upper or lower bound.  If
  all the observations share the same set of common boundaries, these
  are labeled as \sQuote{L1} and \sQuote{U1} for lower and upper
  boundary of the first interval, and so forth.
}
\usage{
boundaries(object)
}
\arguments{
  \item{object}{
    and 'intReg' object
}
}
\value{
  a numeric vector, list of all interval boundaries used in the model in
  increasing order
}
\author{
  Ott Toomet 
}
\keyword{methods}
\keyword{utilities}
