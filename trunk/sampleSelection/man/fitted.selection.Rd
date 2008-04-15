\name{fitted.selection}
\alias{fitted.selection}

\title{Fitted Values of Selection Models}

\description{
   Calculate fitted values of sample selection models
}

\usage{
\method{fitted}{selection}(object, part = "outcome", ... )
}

\arguments{
   \item{object}{object of class \code{selection}.}
   \item{part}{character string indication which fitted values to extract:
      "outcome" for the fitted values of the outcome equation or
      "selection" for the fitted values of the selection equation.}
   \item{\dots}{further arguments passed to other methods
      (e.g. \code{\link{fitted.probit}} or \code{\link{fitted}}).}
}
\value{
   A numeric vector of the fitted values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{selection}},
   \code{\link{fitted}}, and \code{\link{fitted.probit}}}

\keyword{methods}
