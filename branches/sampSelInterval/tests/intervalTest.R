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
hist( dat$yOu )
dat$yO <- cut( dat$yOu, bound )

YS <- dat$yS
XS <- cbind( 1, dat$x1, dat$x2 )
YO <- as.numeric( dat$yO )
XO <- cbind( 1, dat$x1 )

start <- c( betaS, betaO, log( sqrt( sigma ) ), atan( rho ) )
   # the correct starting value of logSigma would be: log( sigma )
names( start ) <- c( "betaS0", "betaS1", "betaS2", "betaO0", "betaO2",
   "logSigma", "atanRho" )

res <- selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound, 
   start = start, printLevel = 1 )

print( res )
print( round( coef( res ), 2 ) )
print( round( coef( summary( res ) ), 2 ) )
print( res$start )
print( summary( res ) )

# add derived coefficients
coefAll <- c( coef( res ),
   sigma = unname( exp( coef( res )[ "logSigma" ] ) ),
   sigmaSq = unname( exp( 2 * coef( res )[ "logSigma" ] ) ),
   rho = unname( tan( coef( res )[ "atanRho" ] ) ) )
print( round( coefAll, 2 ) )

# jacobian
jac <- cbind( diag( length( coef( res ) ) ),
   matrix( 0, length( coef( res ) ), 3 ) )
rownames( jac ) <- names( coef( res ) )
colnames( jac ) <- c( names( coef( res ) ), "sigma", "sigmaSq", "rho" )
jac[ "logSigma", "sigma" ] <- exp( coef( res )[ "logSigma" ] )
jac[ "logSigma", "sigmaSq" ] <- 2 * exp( 2 * coef( res )[ "logSigma" ] )
jac[ "atanRho", "rho" ] <- 1 + ( tan( coef( res )[ "atanRho" ] ) )^2
vcovAll <- t( jac ) %*% vcov( res ) %*% jac
print( round( vcovAll, 2 ) )
print( round( cov2cor( vcovAll ), 2 ) )
print( coefTable( coefAll, sqrt( diag( vcovAll ) ) ) ) 

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
all.equal( as.data.frame( gradStart ), as.data.frame( gradStartNum ) )
library( "miscTools" )
for( i in 1:ncol( gradStart ) ) {
   compPlot( gradStart[ , i ], gradStartNum[ , i ], main = names( start )[i],
      col = ifelse( !YS, "black",
         ifelse( YO == 1, "blue",
            ifelse( YO == 2, "blueviolet",
               ifelse( YO == 3, "red", "green" ) ) ) ) )
}

# log-likelihood values and their gradients at estimated parameters
logLikEst <- intLogLik( coef( res ) )
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
all.equal( resMl[ !names( resMl ) %in% c( "call", "control" ) ],
   resMl2[ !names( resMl2 ) %in% c( "call", "control" ) ] )


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

### tests with Mroz data
data("Mroz87")

## tests with different boundaries
Mroz87$wage_5interval <- cut(Mroz87$wage, br=c(0,2.0,4.0,6.0,8.0,Inf),
   labels=c(1,2,3,4,5))
bound6 <- c(0,2.0,4.0,6.0,8.0,Inf)
Mroz87$wage_6interval <- cut(Mroz87$wage, br=c(0,2.0,4.0,6.0,8.0,10.0,Inf),
   labels=c(1,2,3,4,5,6))
bound7 <- c(0,2.0,4.0,6.0,8.0,10.0,Inf)
Mroz87$wage_7interval <- cut(Mroz87$wage, br=c(0,2.0,4.0,6.0,8.0,10.0,12.0,Inf),
   labels=c(1,2,3,4,5,6,7))
bound8 <- c(0,2.0,4.0,6.0,8.0,10.0,12.0,Inf)

Wage5 <- selection( lfp ~ huswage + kids5 + mtr + fatheduc + educ + city, 
   wage_5interval ~ educ + city, data = Mroz87, boundaries = bound6 )
print(summary(Wage5))

Wage6 <- selection( lfp ~ huswage + kids5 + mtr + fatheduc + educ + city, 
   wage_6interval ~ educ + city, data = Mroz87, boundaries = bound7 )
print(summary(Wage6))

Wage7 <- selection( lfp ~ huswage + kids5 + mtr + fatheduc + educ + city, 
   wage_7interval ~ educ + city, data = Mroz87, boundaries = bound8 )
print(summary(Wage7))

## tests with different specifications
# low number of variables - NA in initial gradient
try(selection( lfp ~ huswage + educ, 
   wage_5interval ~ educ, data = Mroz87, boundaries = bound6))

# adding wife's marginal tax rate (mtr)
spec1 <- selection( lfp ~ huswage + educ + mtr, 
   wage_5interval ~ educ, data = Mroz87, boundaries = bound6)
print(summary(spec1))

# adding continuous variables only
spec2 <- selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, boundaries = bound6)
print(summary(spec2))

# adding dummy variables (city, huscoll)
spec3 <- selection( lfp ~ huswage + mtr + fatheduc + educ + city + huscoll, 
   wage_5interval ~ educ + exper + city, data = Mroz87, boundaries = bound6)
print(summary(spec3))

# only dummy variables as indepdent variables
spec4 <- selection( lfp ~ city + wifecoll, 
   wage_5interval ~ wifecoll, data = Mroz87, boundaries = bound6)
print(summary(spec4))

## trying lrtest and waldtest
library(lmtest)
lrtest(spec1,spec2)
lrtest(spec2,spec3)
waldtest(spec1,spec2,spec3)

## Trying different start value inputs
try(selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, 
   boundaries = bound6, start = c(1,-1,1,-1,-1,-1,1,1,0.5,-0.1) ))
warnings()

try(selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, 
   boundaries = bound6, start = c(3,-0.5,0.2,-4,-0.1,-0.5,0.1,0.2,0.5,-0.1) ))
warnings()

try(selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, 
   boundaries = bound6, start = c(3,-0.3,0.1,-5,-0.1,-0.4,0.2,0.1,0.5,-0.1) ))

## Testing estimations with 'problems'
# NAs in independent variables
Mroz87$huswage[Mroz87$huswage < 4] <- NA 
try(selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, boundaries = bound6))
Mroz87$educ[Mroz87$educ > 12] <- NA
try(selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, boundaries = bound6))

# NAs in dependent variable
Mroz87$wage_5interval[Mroz87$educ < 12] <- NA
try(selection( lfp ~ huswage + educ + mtr + fatheduc, 
   wage_5interval ~ educ + exper, data = Mroz87, boundaries = bound6))


