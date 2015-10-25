## The comprehensive test suite
## Not to be included in CRAN version
## Hence may run long

## Example of observation-specific boundaries
## Estimate the willingness to pay for the Kakadu National Park
## Data given in intervals -- 'lower' for lower bound and 'upper' for upper bound.
## Note that dichotomous-coice answers are already coded to 'lower' and 'upper'
set.seed(1)
options(digits=3)
library(intReg)

## Test a mixed interval complex model and methods
data(Kakadu, package="Ecdat")
## Estimate in log form, change 999 to Inf
Kakadu <- Kakadu[sample(nrow(Kakadu), 100),]
                           # Speed up the tests
lb <- log(Kakadu$lower)
ub <- Kakadu$upper
ub[ub > 998] <- Inf
ub <- log(ub)
## Artifically create a few point observations
iP <- sample(nrow(Kakadu), 15)
iP <- iP[!is.infinite(lb[iP])]
ub[iP] <- lb[iP]
##
y <- cbind(lb, ub)
m <- intReg(y ~ sex + log(income) + age + schooling + 
              recparks + jobs + lowrisk + wildlife + future + aboriginal + finben +
              mineparks + moreparks + gov +
              envcon + vparks + tvenv + major, data=Kakadu)
## You may want to compare the results to Werner (1999),
## Journal of Business and Economics Statistics 17(4), pp 479-486

cat("Coefficients:\n")
print(coef(m))
cat("Coefficients, with boundaries:\n")
print(coef(m, boundaries=TRUE))
cat("stdEr:\n")
print(stdEr(m))
cat("stdEr, with boundaries:\n")
print(stdEr(m, boundaries=TRUE))
cat("Summary:\n")
print(summary(m))
cat("Summary, with boundaries:\n")
print(summary(m, boundaries=TRUE))

## test model.matrix
mm <- model.matrix(m)
cat("Model matrix (sample):\n")
print(mm[1:10,])
## Test model.frame
mf <- model.frame(m)
cat("Model frame (sample):\n")
print(mf[1:10,])
cat("Model response (sample):\n")
print(model.response(mf)[1:10])

## test utility functions
cat("Boundaries:\n")
print(boundaries(m))
cat("Disturbances:\n")
print(disturbances(m))
cat("Intervals (sample):\n")
print(intervals(m)[sample(seq(nObs(m)), 10),])
cat("intervalObs:\n")
print(intervalObs(m))

##
## Example of common intervals for all the observations
##
cat("Common intervals example:\n")
data(Bwages, package="Ecdat")
Bwages <- Bwages[sample(nrow(Bwages), 200),]
## calculate an ordinary Mincer-style wage regression.  
## Note: gross hourly wage rate in EUR
intBound <- c(0, 5, 10, 15, 25, Inf)
salary <- cut(Bwages$wage, intBound)
m <- intReg(salary ~ factor(educ) + poly(exper, 2), data=Bwages,
            boundaries=log(intBound))
## Note: use logs for the intervals in Euros.  We do not have to
## transform salaris to log form as this does not change the intervals.
## Ignore any warnings
cat("Summary, common boundaries:\n")
print(summary(m))

## test utility functions for common intervals
cat("Boundaries:\n")
print(boundaries(m))
cat("Intervals (sample):\n")
print(intervals(m)[1:10,])

## Test model.response
cat("model response (sample):\n")
print(model.response(mf)[1:10])

## test predictions
Ey <- predict(m, type="link")
cat("Link prediction (sample):\n")
print(Ey[1:10])
Eyc <- predict(m, type="linkConditional")
cat("Conditional mean prediction (sample):\n")
print(Eyc[1:10])
## predictions with new data
yInt <- cut(rnorm(5), breaks=c(-4, -3, -2, -1, 0, 1, 2, 3, 4))
newdat <- data.frame(educ=sample(levels(factor(Bwages$educ)), 5),
                     exper=runif(5, 0, 10))
                           # newdat only includes: should work with 'link'
cat("Predicted link:\n")
print(predict(m, newdata=newdat, type="link"))
cat("Predicted expected value conditional on interval:\n")
print(try(predict(m, newdata=newdat, type="linkConditional")))
                           # Error: must include intervals for 'linkConditional'
newdat <- cbind(yInt=yInt, newdat)
print(try(predict(m, newdata=newdat, type="linkConditional")))
                           # should work

##
## Small data, large number of intervals (by Thierry Kalisa)
##
a <- c(0.002300, 0.020000, 0.000150, 0.000005, 0.002300, 0.000045, 0.000150,
       0.000110, 0.000110, 0.000005, 0.010000, 0.000490, 0.000110, 0.000005,
       0.000600, 0.000380, 0.000600, 0.005275, 0.005275, 0.000045, 0.000075,
       0.000600, 0.000600, 0.005275, 0.000075, 0.001650, 0.001100, 0.000005,
       0.000025, 0.005275, 0.000150, 0.005275, 0.000005, 0.000110, 0.000270,
       0.000600, 0.000600, 0.000380, 0.000110, 0.000380, 0.000270, 0.000490,
       0.000045, 0.000110, 0.000110, 0.000150, 0.000005, 0.000110, 0.000045,
       0.005275, 0.000600, 0.000200, 0.003475, 0.005275, 0.000005, 0.000600,
       0.000200, 0.000075, 0.000600, 0.000600, 0.000075, 0.000230, 0.000490,
       0.005275, 0.000230, 0.000110, 0.000490, 0.000045, 0.000075, 0.001650,
       0.000600, 0.000490, 0.000005, 0.003475, 0.001650, 0.000150, 0.000380,
       0.017500, 0.003475, 0.000270, 0.000230, 0.005275, 0.000045, 0.000045,
       0.000075, 0.003475, 0.000150, 0.002300, 0.001650, 0.001100, 0.000005,
       0.000075, 0.000025, 0.000025, 0.000150, 0.001100)
b <- c(0.003475, 0.040000, 0.005275, 0.040000, 0.015000, 0.001100, 0.000380,
       0.003475, 0.003475, 0.040000, 0.020000, 0.007075, 0.000490, 0.003475,
       0.007075, 0.005275, 0.012500, 0.012500, 0.010000, 0.000270, 0.000200,
       0.002300, 0.010000, 0.010000, 0.001650, 0.003475, 0.005275, 0.003475,
       0.003475, 0.010000, 0.000600, 0.020000, 0.000045, 0.001650, 0.010000,
       0.005275, 0.020000, 0.001650, 0.005275, 0.003475, 0.003475, 0.007075,
       0.002300, 0.010000, 0.000270, 0.000270, 0.003475, 0.000600, 0.000270,
       0.007075, 0.003475, 0.010000, 0.010000, 0.012500, 0.000045, 0.010000,
       0.003475, 0.010000, 0.012500, 0.003475, 0.000380, 0.003475, 0.005275,
       0.008650, 0.000600, 0.002300, 0.003475, 0.005275, 0.003475, 0.003475,
       0.003475, 0.002300, 0.000025, 0.017500, 0.005275, 0.003475, 0.001650,
       0.020000, 0.040000, 0.001650, 0.003475, 0.008650, 0.000200, 0.000110,
       0.000490, 0.040000, 0.000600, 0.020000, 0.005275, 0.008650, 0.000490,
       0.005275, 0.000230, 0.000200, 0.000270, 0.005275)
c <-c(3, 4, 3, 3, 3, 1, 2, 1, 3, 4, 2, 2, 1, 2, 1, 2, 2, 1, 3, 2, 2, 3, 1, 2, 1, 2, 3, 2, 4, 3, 4, 2,
      4, 2, 1, 2, 4, 3, 2, 3, 2, 2, 3, 4, 2, 1, 3, 3, 1, 1, 2, 1, 2, 2, 1, 3, 1, 1, 2, 3, 2, 2, 3, 1,
      3, 2, 2, 1, 2, 2, 2, 2, 1, 3, 2, 3, 2, 1, 1, 2, 2, 1, 1, 2, 3,
      1, 2, 3, 2, 2, 1, 1, 4, 1, 3, 3)
ab <- cbind(a,b)
mNorm <- intReg(ab~c)
print(summary(mNorm))

## Test the same with cloglog disturbances
m <- intReg(ab~c, method="cloglog")
print(summary(m))
cat("Disturbances:\n")
print(disturbances(m))

## Test precision of intervals
