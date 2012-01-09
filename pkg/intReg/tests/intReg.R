## Example of observation-specific boundaries
## Estimate the willingness to pay for the Kakadu National Park
## Data given in intervals -- 'lower' for lower bound and 'upper' for upper bound.
## Note that dichotomous-coice answers are already coded to 'lower' and 'upper'
set.seed(1)
options(digits=5)
library(intReg)
library(Ecdat)
data(Kakadu)
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
## You may want to compare the results to Werner (1999),
## Journal of Business and Economics Statistics 17(4), pp 479-486
print(summary(m))

## test model.matrix
mm <- model.matrix(m)
print(mm[i <- sample(nrow(mm), 10),])

##
## Example of common intervals for all the observations
##
library(Ecdat)
data(Bwages)
## calculate an ordinary Mincer-style wage regression.  
## Note: gross hourly wage rate in EUR
intervals <- c(0, 5, 10, 15, 25, Inf)
salary <- cut(Bwages$wage, intervals)
int <- intReg(salary ~ factor(educ) + poly(exper, 2), data=Bwages, boundaries=log(intervals))
## Note: use logs for the intervals in Euros.  We do not have to
## transform salaris to log form as this does not change the intervals.
## Ignore any warnings
cat("Interval regression:\n")
print(summary(int))

## Test model.frame
mf <- model.frame(int)
print(mf[i <- sample(nrow(mf), 10),])
print(model.response(mf)[i])
