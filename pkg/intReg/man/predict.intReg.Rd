\name{predict.intReg}
\alias{predict.intReg}
\title{
  predicted outcomes of 'intReg' models
}
\description{
  predict the latent link value for interval regression models, either
  unconditional or conditional of being in the observed interval.
}
\usage{
\method{predict}{intReg}(object, newdata = NULL, type = "link", ...)
}
\arguments{
  \item{object}{
    \dQuote{intReg} object, estimated interval regression model
}
  \item{newdata}{
    data frame for new data
}
  \item{type}{
    what to predict: \dQuote{link} for the expected
    latent link value,
    \dQuote{linkConditional} the expected latent value conditional of
    being in the observed interval
}
\item{\dots}{
  other arguments for the model
}
}
\value{
  a numeric vector, the predicted values by observation
}
\author{
  Ott Toomet
}

\seealso{
\code{\link{intReg}}
}
\examples{
}
\keyword{methods}
\keyword{models}
