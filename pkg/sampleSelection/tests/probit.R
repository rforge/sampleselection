library( "sampleSelection" )
options( digits = 3 )

## loading and preparing data
data( Mroz87 )
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
Mroz87$age30.39 <- Mroz87$age < 40
Mroz87$age50.60 <- Mroz87$age >= 50

## A simple single MC trial: note probit assumes normal errors
set.seed( 20080225 )
x <- runif( 100 )
e <- rnorm( 100 )
y <- 2 * x + e
probitResult <- probit( (y > 0) ~ x )
print( probitResult )
summary( probitResult )
nObs( probitResult )
df.residual( probitResult )
logLik( probitResult )
fitted( probitResult )
residuals( probitResult, type = "response" )
residuals( probitResult, type = "pearson" )
residuals( probitResult, type = "deviance" )
all.equal( residuals( probitResult, type = "deviance" ),
   residuals( probitResult ) )

# estimation with glm()
probitResult2 <- glm( (y > 0) ~ x, family = binomial( link = "probit" ) )
all.equal( coef( probitResult ), coef( probitResult2 ), tol = 1e-4 )
all.equal( stdEr( probitResult ), stdEr( probitResult2 ), tol = 1e-1 )
all.equal( logLik( probitResult ), logLik( probitResult2 ) )
all.equal( fitted( probitResult ), fitted( probitResult2 ), tol = 1e-4 )
all.equal( residuals( probitResult, type = "response" ),
   residuals( probitResult2, type = "response" ), tol = 1e-4 )
all.equal( residuals( probitResult, type = "pearson" ),
   residuals( probitResult2, type = "pearson" ), tol = 1e-4 )
all.equal( residuals( probitResult, type = "deviance" ),
   residuals( probitResult2, type = "deviance" ), tol = 1e-4 )
all.equal( residuals( probitResult ), residuals( probitResult2 ), tol = 1e-4 )

# estimation with equal weights
we <- rep( 0.5, 100 )
probitResultWe <- probit( (y > 0) ~ x, weights = we )
summary( probitResultWe )
nObs( probitResultWe )
df.residual( probitResultWe )
logLik( probitResult )
all.equal( coef( probitResult ), coef( probitResultWe ), tol = 1e-4 )
all.equal( logLik( probitResult ) * 0.5, logLik( probitResultWe ) * 1 )

# estimation with equal weights with glm()
probitResultWe2 <- glm( (y > 0) ~ x, family = binomial( link = "probit" ),
   weights = we )
all.equal( coef( probitResultWe ), coef( probitResultWe2 ), tol = 1e-4 )
all.equal( stdEr( probitResultWe ), stdEr( probitResultWe2 ), tol = 1e-1 )
logLik( probitResultWe2 )
all.equal( fitted( probitResultWe ), fitted( probitResultWe2 ), tol = 1e-4 )
all.equal( residuals( probitResultWe, type = "response" ),
   residuals( probitResultWe2, type = "response" ), tol = 1e-4 )
all.equal( residuals( probitResultWe, type = "pearson" ),
   residuals( probitResultWe2, type = "pearson" ), tol = 1e-4 )
all.equal( residuals( probitResultWe, type = "deviance" ),
   residuals( probitResultWe2, type = "deviance" ), tol = 1e-4 )

# estimation with weights to account for stratified sampling
# proportion in the "population"
yProbPop <- sum( y > 0 ) / length( y )
yProbPop
# stratified sample with all observations with y = 0
sampStrat <- y <= 0 | rnorm( length( y ) ) > 0.25
sum( sampStrat )
# stratified sample of y and x
ySamp <- y[ sampStrat ]
xSamp <- x[ sampStrat ]
# proportion in the "sample"
yProbSamp <- sum( ySamp > 0 ) / length( ySamp )
yProbSamp
# unweighted estimation (ignoring the stratification)
probitResultStrat <- probit( (ySamp > 0) ~ xSamp )
summary( probitResultStrat )
nObs( probitResultStrat )
df.residual( probitResultStrat )
logLik( probitResultStrat )
# weights
wStrat <- ifelse( ySamp > 0, yProbPop / yProbSamp,
   ( 1 - yProbPop ) / ( 1 - yProbSamp ) )
probitResultStratW <- probit( (ySamp > 0) ~ xSamp, weights = wStrat )
summary( probitResultStratW )
nObs( probitResultStratW )
df.residual( probitResultStratW )
logLik( probitResultStratW )

# estimation with weights to account for stratified sampling with glm()
probitResultStratW2 <- glm( (ySamp > 0) ~ xSamp, 
   family = binomial( link = "probit" ), weights = wStrat )
all.equal( coef( probitResultStratW ), coef( probitResultStratW2 ), tol = 1e-4 )
all.equal( stdEr( probitResultStratW ), stdEr( probitResultStratW2 ), tol = 1e-1 )
logLik( probitResultStratW2 )
all.equal( fitted( probitResultStratW ), fitted( probitResultStratW2 ), tol = 1e-4 )
all.equal( residuals( probitResultStratW, type = "response" ),
   residuals( probitResultStratW2, type = "response" ), tol = 1e-4 )
all.equal( residuals( probitResultStratW, type = "pearson" ),
   residuals( probitResultStratW2, type = "pearson" ), tol = 1e-4 )
all.equal( residuals( probitResultStratW, type = "deviance" ),
   residuals( probitResultStratW2, type = "deviance" ), tol = 1e-4 )


## female labour force participation probability
lfpResult <- probit( lfp ~ kids + age30.39 + age50.60 + educ + hushrs +
   huseduc + huswage + mtr + motheduc, data = Mroz87 )
print( lfpResult )
summary( lfpResult )
nObs( lfpResult )
df.residual( lfpResult )
logLik( lfpResult )
fitted( lfpResult )
residuals( lfpResult, type = "response" )
residuals( lfpResult, type = "pearson" )
residuals( lfpResult, type = "deviance" )

# estimation with glm()
lfpResult2 <- glm( lfp ~ kids + age30.39 + age50.60 + educ + hushrs +
      huseduc + huswage + mtr + motheduc, data = Mroz87, 
   family = binomial( link = "probit" ) )
all.equal( coef( lfpResult ), coef( lfpResult2 ), tol = 1e-3 )
all.equal( stdEr( lfpResult ), stdEr( lfpResult2 ), tol = 1e-1 )
all.equal( logLik( lfpResult ), logLik( lfpResult2 ) )
all.equal( fitted( lfpResult ), fitted( lfpResult2 ), tol = 1e-4 )
all.equal( residuals( lfpResult, type = "response" ),
   residuals( lfpResult2, type = "response" ), tol = 1e-4 )
all.equal( residuals( lfpResult, type = "pearson" ),
   residuals( lfpResult2, type = "pearson" ), tol = 1e-4 )
all.equal( residuals( lfpResult, type = "deviance" ),
   residuals( lfpResult2, type = "deviance" ), tol = 1e-4 )
all.equal( residuals( lfpResult ), residuals( lfpResult2 ), tol = 1e-4 )


## Greene( 2003 ): example 22.8, page 786 (only probit part )
greene <- probit( lfp ~ age + I( age^2 ) + faminc + kids + educ, data = Mroz87 )
print( greene )
summary( greene )
nObs( greene )
df.residual( greene )
logLik( greene )
fitted( greene )
residuals( greene, type = "response" )
residuals( greene, type = "pearson" )
residuals( greene, type = "deviance" )

# estimation with glm()
greene2 <- glm( lfp ~ age + I( age^2 ) + faminc + kids + educ, 
   data = Mroz87, family = binomial( link = "probit" ) )
all.equal( coef( greene ), coef( greene2 ), tol = 1e-3 )
all.equal( stdEr( greene ), stdEr( greene2 ), tol = 1e-1 )
all.equal( logLik( greene ), logLik( greene2 ) )
all.equal( fitted( greene ), fitted( greene2 ), tol = 1e-4 )
all.equal( residuals( greene, type = "response" ),
   residuals( greene2, type = "response" ), tol = 1e-4 )
all.equal( residuals( greene, type = "pearson" ),
   residuals( greene2, type = "pearson" ), tol = 1e-4 )
all.equal( residuals( greene, type = "deviance" ),
   residuals( greene2, type = "deviance" ), tol = 1e-4 )
all.equal( residuals( greene ), residuals( greene2 ), tol = 1e-4 )


## factors as dependent variable (from Achim Zeileis)
probit( lfp ~ exper, data = Mroz87 )
probit( factor( lfp ) ~ exper, data = Mroz87 )
probit( factor( lfp, labels = c( "no", "yes" ) ) ~ exper, data = Mroz87 )

## NA in data/linear predictors/na.exclude (by Gabor Grothendieck)
x <- runif(20)
y <- x + rnorm(length(x)) > 0
y[1] <- y[4] <- NA
m <- probit(y ~ x, na.action = na.exclude)
length(linearPredictors(m))
                           # should be 20
nObs( m )
df.residual( m )
logLik( m )

# example, where the "cutoff" in the log likelihood function is used
# (inspired by an example from Jon K. Peck)
set.seed( 123 )
nObs <- 1000
dta2 <- data.frame( id = c( 1:nObs ) )
for( i in 1:3 ) {
   dta2[[ paste( "x", i, sep = "" ) ]] <- rnorm( nObs, 5, 3 )
}
dta2$y <- with( dta2, -5.5 + 0.25 * x1 + 0.6 * x2 + 0.85 * x3 + rnorm( nObs ) ) > 0
p2 <- probit( y ~ x1 + x2 + x3, data=dta2 )
nObs( p2 )
df.residual( p2 )
logLik( p2 )

## This test probes the extreme outliers.  We generate a normal
## model, and add a positive and a negative extreme outlier.  The
## likelihood code should be robust and not crash.  The estimates
## (and standard errors) are probably of little value.
x <- runif(20) - 0.5
y <- x + rnorm(20) > 0
x[1] <- 1000
y[1] <- FALSE
x[2] <- -1000
y[2] <- TRUE
m <- probit(y ~ x)
print(summary(m))
nObs( m )
df.residual( m )
logLik( m )
