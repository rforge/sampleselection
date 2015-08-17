\name{boundaries}
\alias{boundaries}
\title{
  Return interval boundaries for intReg models 
}
\description{
  Return a vector of all interval boundaries in increasing order.  This
  also includes eventual \code{-Inf} and \code{Inf} if used.  No
  distinction is made if a boundary is used as upper or lower bound.  In
  case of observation-specific intervals, the list may be long.
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
