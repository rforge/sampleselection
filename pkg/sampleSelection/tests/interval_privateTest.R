library( "sampleSelection" )
library( "mvtnorm" )
options( digits = 2 )

nObs <- 300
betaS <- c( 1, 1, -1 )
betaO <- c( 10, 4 )
rho <- 0.4
sigma <- 5
bound <- c(-Inf,5,15,Inf)


set.seed(123)
dat <- data.frame( x1 = rnorm( nObs ), x2 = rnorm( nObs ) )
vcovMat <- matrix( c( 1, rho*sigma, rho*sigma, sigma^2 ), nrow = 2 )
eps <- rmvnorm( nObs, sigma = vcovMat )
dat$epsS <- eps[,1]
dat$epsO <- eps[,2]
dat$yS <- with( dat, betaS[1] + betaS[2] * x1 + betaS[3] * x2 + epsS ) > 0
# table(dat$yS)
dat$yOu <- with( dat, betaO[1] + betaO[2] * x1 + epsO )
dat$yOu[ !dat$yS ] <- NA
# hist( dat$yOu )
dat$yO <- cut( dat$yOu, bound )

YS <- dat$yS
XS <- cbind( 1, dat$x1, dat$x2 )
YO <- as.numeric( dat$yO )
XO <- cbind( 1, dat$x1 )


start <- c( betaS, betaO, log( sqrt( sigma ) ), atanh( rho ) )
# the correct starting value of logSigma would be: log( sigma )
names( start ) <- c( "betaS0", "betaS1", "betaS2", "betaO0", "betaO2",
   "logSigma", "atanhRho" )

res <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = start, printLevel = 1 )

print( res )
print( round( coef( res ), 2 ) )
print( round( coef( summary( res ) ), 2 ) )
print( res$start )
print( summary( res ) )

# with derived coefficients
print( round( res$coefAll, 2 ) )
print( round( res$vcovAll, 2 ) )
print( round( cov2cor( res$vcovAll ), 2 ) )
print( coefTable( res$coefAll, sqrt( diag( res$vcovAll ) ) ) ) 

maxLik:::summary.maxLik( res )

# function that returns log-likelihood values
intLogLik <- function( param ) {
   ll <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
      start = param, returnLogLikStart = TRUE )
   return( ll )
}


# log-likelihood values and their gradients at starting values of parameters
logLikStart <- intLogLik( start )
print( c( logLikStart ) )
gradStart <- attr( logLikStart, "gradient" )
colnames( gradStart ) <- names( start )
print( gradStart )
# numeric gradients
gradStartNum <- numericGradient( intLogLik, t0 = start )
colnames( gradStartNum ) <- names( start )
all.equal( as.data.frame( gradStart ), as.data.frame( gradStartNum ),
   tol = 1e-5 )
# library( "miscTools" )
# for( i in 1:ncol( gradStart ) ) {
#    compPlot( gradStart[ , i ], gradStartNum[ , i ], main = names( start )[i],
#       col = ifelse( !YS, "black",
#          ifelse( YO == 1, "blue",
#             ifelse( YO == 2, "blueviolet",
#                ifelse( YO == 3, "red", "green" ) ) ) ) )
# }

# log-likelihood values and their gradients at estimated parameters
logLikEst <- intLogLik( coef( res, part = "est" ) )
print( c( logLikEst ) )
gradEst <- gradStart <- attr( logLikEst, "gradient" )
colnames( gradEst ) <- names( start )
print( gradEst )


# tests with automatically generated starting values (ML estimation)
resMl <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = "ml", printLevel = 1 )
print( resMl )
print( round( coef( resMl ), 2 ) )
print( round( coef( summary( resMl ) ), 2 ) )
print( resMl$start )
print( summary( resMl ) )

resMl2 <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound )
all.equal( resMl[ !names( resMl ) %in% c( "call", "control", "objectiveFn" ) ],
   resMl2[ !names( resMl2 ) %in% c( "call", "control", "objectiveFn" ) ] )


# tests with automatically generated starting values (2-step estimation)
res2s <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = "2step", printLevel = 1 )
print( res2s )
print( round( coef( res2s ), 2 ) )
print( round( coef( summary( res2s ) ), 2 ) )
print( res2s$start )
print( summary( res2s ) )


# tests with incorrectly specified starting values
try( selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = "wrong" ) )
try( selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = rep( 1, 11 ) ) )

# tests with incorrectly specified boundaries
try( selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = 1:6, 
   start = start ) )
try( selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = 4:1, 
   start = start ) )


# Test estimation with empty interval and yO as integer
bound <- c(-Inf,4.9,5,15,Inf)
dat$yO <- cut( dat$yOu, br = bound )
empty1 <- selection( yS ~ x1 + x2, as.integer(yO) ~ x1, data = dat,
   boundaries = bound, start = start, printLevel = 1 )

# Test estimation with empty interval and yO as numeric variable
empty2 <- selection( yS ~ x1 + x2, as.numeric(yO) ~ x1, data = dat,
   boundaries = bound, start = start, printLevel = 1 )
all.equal( coef(empty1), coef(empty2) )

# Test estimation with empty interval and yO as factor
empty3 <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = start, printLevel = 1 )
all.equal( coef(empty1), coef(empty3) )


## Testing estimations with NAs
# NAs in independent variables
dat$x1[dat$x1 > 0.1 & dat$x1 < 1.2] <- NA
dat$x2[dat$x2 < -0.5 & dat$x2 > -0.7] <- NA
NAres1 <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = start, printLevel = 1 )
print(NAres1)

# NAs in dependent variable
dat$yS[dat$x1 > 1.1 & dat$x1 < 1.4] <- NA
NAres2 <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = start, printLevel = 1 )
print(NAres2)


### tests with Mroz data
data("Mroz87")
bounds <- c(0,2.0,4.0,6.0,8.0,Inf)
Mroz87$wage_5interval <- cut( Mroz87$wage, br = bounds )

## tests with different specifications
# low number of variables
spec1 <- selection( lfp ~ huswage + educ, 
   wage_5interval ~ huswage, data = Mroz87, boundaries = bounds)
print(summary(spec1))

# adding continuous variables only
spec2 <- selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, boundaries = bounds)
print(summary(spec2))

# adding dummy variables (city, huscoll)
spec3 <- selection( lfp ~ huswage + mtr + fatheduc + educ + city + huscoll, 
   wage_5interval ~ educ + exper + city, data = Mroz87, boundaries = bounds)
print(summary(spec3))

# only dummy variables as indepdent variables
spec4 <- selection( lfp ~ city + wifecoll, 
   wage_5interval ~ city, data = Mroz87, boundaries = bounds)
print(summary(spec4))

# trying lrtest and waldtest
library(lmtest)
lrtest(spec1,spec2)
lrtest(spec2,spec3)
waldtest(spec1,spec2,spec3)


### tests with Smoke data
data(Smoke)

## tests with different specifications
bounds <- c(0,5,10,20,50,Inf)

# test with low number of variables
Smoke_spec1 <- selection( smoker ~ educ + age, 
   cigs_intervals ~ educ, data = Smoke, boundaries = bounds)

# test with more variables
Smoke_spec2 <- selection( smoker ~ educ + age + restaurn, 
   cigs_intervals ~ educ + income + restaurn, data = Smoke, boundaries = bounds)

# testing models against each other
lrtest(Smoke_spec1, Smoke_spec2)
waldtest(Smoke_spec1, Smoke_spec2)