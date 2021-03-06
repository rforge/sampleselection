library( sampleSelection )
options( digits = 6 )

## loading and preparing data
data( Mroz87 )
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
Mroz87$age30.39 <- Mroz87$age < 40
Mroz87$age50.60 <- Mroz87$age >= 50

## A simple single MC trial: note probit assumes normal errors
set.seed( 20080225 )
x <- runif( 100 )
e <- 0.5 * rnorm( 100 )
y <- x + e
probitResult <- probit( (y > 0) ~ x )
print( probitResult )
summary( probitResult )
fitted( probitResult )
residuals( probitResult, type = "response" )
residuals( probitResult, type = "pearson" )
residuals( probitResult, type = "deviance" )

## female labour force participation probability
lfpResult <- probit( lfp ~ kids + age30.39 + age50.60 + educ + hushrs +
   huseduc + huswage + mtr + motheduc, data = Mroz87 )
print( lfpResult )
summary( lfpResult )
fitted( lfpResult )
residuals( lfpResult, type = "response" )
residuals( lfpResult, type = "pearson" )
residuals( lfpResult, type = "deviance" )

## Greene( 2003 ): example 22.8, page 786 (only probit part )
greene <- probit( lfp ~ age + I( age^2 ) + faminc + kids + educ, data = Mroz87 )
print( greene )
summary( greene )
fitted( greene )
residuals( greene, type = "response" )
residuals( greene, type = "pearson" )
residuals( greene, type = "deviance" )

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
