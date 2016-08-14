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
    a list or frame for new data.  If \code{type} is
    \sQuote{linkConditional}, the new data frame must include new
    intervals.  These can either be supplied through two variables
    \sQuote{LB} and \sQuote{UB} for lower and upper bound, or a single
    variable \sQuote{yInt} which is a factor with
    \code{\link{cut}}--style interval coding in the form
    \dQuote{[LB,UB)}.  (The brackets and parentheses are equivalent.)
  }
  \item{type}{
    character, what to predict:
    \itemize{
      \item{\dQuote{link}} \eqn{y^*}{y*} for the expected
      latent link value,
      \item{\dQuote{linkConditional}} the expected latent value conditional of
      being in the observed interval \eqn{(l,u)}, \eqn{\mathrm{E} [y^*| l \le y^* < u]}{E[y*|l
	<= y* < u]}
    }
  }
  \item{\dots}{
    currently ignored
  }
}
\value{
  a numeric vector, the predicted values
}
\author{
  Ott Toomet
}

\seealso{
\code{\link{intReg}}
}
\examples{
   data(Kakadu, package="Ecdat")
   set.seed(1)
   Kakadu <- Kakadu[sample(nrow(Kakadu), 400),]
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
   ## New example with new data.  Note the intervals are only necessary for
   ## 'linkConditional' prediction
   ## For simplicity, sample from the original one
   newDat <- Kakadu[sample(nrow(Kakadu), 10),]
   ## There are two ways to create new intervals:
   ## either varaibles 'LB' and 'UB' in the new data frame:
   newdatA <- newDat
   newdatA$LB <- log(newDat$lower)
   newdatA$UB <- log(newDat$upper)
   EycA <- predict(m, newdata=newdatA, type="linkConditional")
   print(head(EycA))
   ## ... or by introducing a variable 'yInt' which is a factor
   ## with 'cut'-style intervals like '[LB,UB)':
   newdatB <- newDat
   y <- rnorm(nrow(newDat), 6, 3)
   newdatB$yInt <- cut(y, breaks=c(0, 2, 5, 20, 50, 100, 250, Inf))
   EycB <- predict(m, newdata=newdatB, type="linkConditional")
   print(head(EycB))
}
\keyword{methods}
\keyword{models}
