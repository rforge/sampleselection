\name{tobitBayes}
\alias{tobitBayes}
\title{Bayesian Tobit Estimation}

\description{
   Estimate a standard Tobit model
   (endogenous variable is left-censored at zero)
   using Bayesian Markov Chain Monte Carlo (MCMC) Gibbs sampling.
   Currently, this is implemented for panel data only.
}

\usage{
tobitBayes( formula, data, nRep = 500, nBurn = 100, ... )
}

\arguments{
   \item{formula}{a symbolic description for the model to be estimated.}
   \item{data}{a panel data frame created with \code{\link[plm]{plm.data}}.}
   \item{nRep}{number of replications (without burn-in replications).}
   \item{nBurn}{number of replications in burn-in period.}
   \item{\ldots}{further arguments passed to \code{\link[plm]{plm}}
      (see details).}
}

\details{
   \code{\link[plm]{plm}} is called with arguments
   \code{model = "random"} and \code{effect = "individual"}.
   This might become more flexible in the future.
}

\value{
   a list of class \code{tobitBayes} containing following elements:
   \item{call}{the matched call.}
}

\references{
Bruno, Giuseppe (2004)
\emph{Limited Dependent Panel Data Models: A Comparative Analysis
of Classical and Bayesian Inference among Econometric Packages},
Computing in Economics and Finance 41,
Society for Computational Economics.
}

\author{Arne Henningsen}

\examples{
   data( nlswork )
   # add information on panel structure; using a subset to save time
   library( plm )
   nw <- plm.data( nlswork[1:498,], c( "idcode", "year" ) )
   estResult <- tobitBayes( wks_ue ~ union + age + grade + not_smsa + south +
      occ_code, data = nw, nBurn = 10, nRep = 50 )
   estResult$coef
}

\keyword{models}
