\name{residuals.probit}
\alias{residuals.probit}

\title{Residuals of probit models}
\description{
   Calculate residuals of \code{\link{probit}} models.
}

\usage{
   \method{residuals}{probit}( object, type = "response", ... )
}

\arguments{
   \item{object}{an object of class \code{probit}.}
   \item{type}{the type of residuals which should be returned.
      The alternatives are: "response" (default), "deviance", and "pearson"
      (see \code{\link[stats]{residuals.glm}}
      and Davison & Snell 1991).}
   \item{\dots}{further arguments (currently ignored).}
}

\value{
   A numeric vector of the residuals.
}

\references{
   Davison, A. C. and Snell, E. J. (1991)
   \emph{Residuals and diagnostics.}
   In: Statistical Theory and Modelling. In Honour of Sir David Cox,
   edited by Hinkley, D. V., Reid, N. and Snell, E. J.,
   Chapman & Hall, London.
}

\author{
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}
}

\seealso{\code{\link[stats]{residuals}}, \code{\link{probit}}.}
\keyword{methods}
