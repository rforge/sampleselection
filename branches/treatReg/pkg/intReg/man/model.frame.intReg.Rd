\name{model.frame.intReg}
\alias{model.frame.intReg}

\title{Data used for estimating interval regression model}

\description{
   Return the variables used for estimating an interval regression model
}

\usage{
\method{model.frame}{intReg}( formula, ... )
}

\arguments{
   \item{formula}{object of class \code{intReg}.}
   \item{\dots}{further arguments passed to other methods
      (e.g. \code{\link{model.frame}}).}
}
\value{
  A data.frame containing all variables used for the estimation.

  The response variable is given as a factor (in case of common
  boundaries) or as two numeric vectors, corresponding to the lower- and
  upper boundary for each observation.
}

\author{Arne Henningsen,
   Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{intReg}}, and \code{\link{model.frame}}}
\keyword{methods}
