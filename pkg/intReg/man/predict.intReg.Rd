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
   library(Ecdat)
   data(Kakadu)
   set.seed(1)
   Kakadu <- Kakadu[sample(nrow(Kakadu), 500),]
                           # subsample to speed up the estimation
   ## Estimate in log form, change 999 to Inf
   lb <- log(Kakadu$lower)
   ub <- Kakadu$upper
   ub[ub > 998] <- Inf
   ub <- log(ub)
   y <- cbind(lb, ub)
   m <- intReg(y ~ sex + log(income) + age + schooling + 
                   recparks + jobs + lowrisk + wildlife + future + aboriginal + finben +
                       mineparks + moreparks + gov +
                           envcon + vparks + tvenv + major, data=Kakadu)
   ## Expect the conditional and unconditional predictions to deviate by
   Ey <- predict(m, type="link")
   print(head(Ey))
   Eyc <- predict(m, type="linkConditional")
   print(head(Eyc))
}
\keyword{methods}
\keyword{models}
