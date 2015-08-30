### A general test suite intended to be added to CRAN
### Test general behavior, do not go into depth.
### 
### Note: there are far more thorough tests on the intReg R-forge page.

## Observation-specific boundaries
## Estimate the willingness to pay for the Kakadu National Park
## Data given in intervals -- 'lower' for lower bound and 'upper' for upper bound.
## Note that dichotomous-coice answers are already coded to 'lower' and 'upper'
set.seed(1)
options(digits=4)
library(intReg)
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
y <- cbind(lb, ub)
m <- intReg(y ~ sex + log(income),
            data=Kakadu)
## You may want to compare the results to Werner (1999),
## Journal of Business and Economics Statistics 17(4), pp 479-486

## Test coef, stdEr, summary with and without boundaries
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
cat("Model matrix (sample):\n")
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
library(Ecdat)
data(Bwages)
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

## test predictions
Ey <- predict(m, type="link")
cat("Link prediction (sample):\n")
print(Ey[1:10])
Eyc <- predict(m, type="linkConditional")
cat("Conditional mean prediction (sample):\n")
print(Eyc[1:10])

## Test the same with cloglog disturbances
m <- intReg(salary ~ factor(educ) + poly(exper, 2), data=Bwages,
            boundaries=log(intBound),
            method="logistic")
cat("Same, logit disturbances:\n")
print(summary(m))
cat("Disturbances:\n")
print(disturbances(m))
