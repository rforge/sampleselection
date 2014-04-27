library( "sampleSelection" )
data( "Mroz87" )
options( digits = 3 )

## Greene( 2003 ): example 22.8, page 786
## 2-step estimation
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
greene <- heckit( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
print( greene )
print(summary( greene ))
print(summary( greene$lm ) )
print( summary( greene$probit ) )
print( greene$sigma )
print( greene$rho )
print( coef( greene ) )
print( coef( greene, part = "outcome" ) )
print( coef( summary( greene ) ) )
print( coef( summary( greene ), part = "outcome" ) )
print( vcov( greene ) )
print( vcov( greene, part = "outcome" ) )
nobs( greene )
nObs( greene )
fitted( greene, part = "outcome" )
fitted( greene, part = "selection" )
residuals( greene, part = "outcome" )
residuals( greene, part = "selection", type = "response" )
print( model.matrix( greene, part = "outcome" ) )
print( model.matrix( greene, part = "selection" ) )
print( model.frame( greene ) )

# ML estimation
greeneMl <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
print( greeneMl )
print(summary( greeneMl ))
print( coef( greeneMl ) )
print( coef( greeneMl, part = "outcome" ) )
print( coef( summary( greeneMl ) ) )
print( coef( summary( greeneMl ), part = "outcome" ) )
print( vcov( greeneMl ) )
print( vcov( greeneMl, part = "outcome" ) )
all.equal( nobs( greene ), nobs( greeneMl ) )
all.equal( nObs( greene ), nObs( greeneMl ) )
fitted( greeneMl, part = "outcome" )
fitted( greeneMl, part = "selection" )
residuals( greeneMl, part = "outcome" )
residuals( greeneMl, part = "selection", type = "response" )
all.equal( model.matrix( greene, part = "outcome" )[,-6],
   model.matrix( greeneMl, part = "outcome" ),
   check.attributes = FALSE )
model.matrix( greeneMl, part = "outcome" )
all.equal( model.matrix( greene, part = "selection" ),
   model.matrix( greeneMl, part = "selection" ),
   check.attributes = FALSE )
model.matrix( greeneMl, part = "selection" )
all.equal( model.frame( greene )[,-c(2,12)], model.frame( greeneMl ),
   check.attributes = FALSE )
model.frame( greeneMl )


## Wooldridge( 2003 ): example 17.5, page 590
## 2-step estimation
data( "Mroz87" )
wooldridge <- heckit( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
print( wooldridge )
print( summary( wooldridge ) )
print( summary( wooldridge$lm ) )
print( summary( wooldridge$probit ) )
print( wooldridge$sigma )
print( wooldridge$rho )
print( coef( wooldridge ) )
print( coef( wooldridge, part = "outcome" ) )
print( coef( summary( wooldridge ) ) )
print( coef( summary( wooldridge ), part = "outcome" ) )
print( vcov( wooldridge ) )
print( vcov( wooldridge, part = "outcome" ) )
nobs( wooldridge )
nObs( wooldridge )
round( fitted( wooldridge, part = "outcome" ), digits = 3 )
round( fitted( wooldridge, part = "selection" ), digits = 3 )
round( residuals( wooldridge, part = "outcome" ), digits = 3 )
round( residuals( wooldridge, part = "selection", type = "response" ), digits = 3 )
print( model.matrix( wooldridge, part = "outcome" ) )
print( model.matrix( wooldridge, part = "selection" ) )
print( model.frame( wooldridge ) )

# ML estimation
wooldridgeMl <- selection( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
      kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
print( wooldridgeMl )
print( summary( wooldridgeMl ) )
print( coef( wooldridgeMl ) )
print( coef( wooldridgeMl, part = "outcome" ) )
print( coef( summary( wooldridgeMl ) ) )
print( coef( summary( wooldridgeMl ), part = "outcome" ) )
print( vcov( wooldridgeMl ) )
print( vcov( wooldridgeMl, part = "outcome" ) )
nobs( wooldridgeMl )
nObs( wooldridgeMl )
round( fitted( wooldridgeMl, part = "outcome" ), digits = 3 )
round( fitted( wooldridgeMl, part = "selection" ), digits = 3 )
round( residuals( wooldridgeMl, part = "outcome" ), digits = 3 )
round( residuals( wooldridgeMl, part = "selection", type = "response" ), digits = 3 )
all.equal( model.matrix( wooldridge, part = "outcome" )[,-5],
   model.matrix( wooldridgeMl, part = "outcome" ),
   check.attributes = FALSE )
model.matrix( wooldridgeMl, part = "outcome" )
all.equal( model.matrix( wooldridge, part = "selection" ),
   model.matrix( wooldridgeMl, part = "selection" ),
   check.attributes = FALSE )
model.matrix( wooldridgeMl, part = "selection" )
all.equal( model.frame( wooldridge )[,-c(2,11)], model.frame( wooldridgeMl ),
   check.attributes = FALSE )
print( model.frame( wooldridgeMl ) )


## Tobit 5 Example from the selection paper
library( "mvtnorm" )
set.seed(0)
vc <- diag(3)
vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
vc[upper.tri(vc)] <- vc[lower.tri(vc)]
eps <- rmvnorm(500, rep(0, 3), vc )
xs <- runif(500)
ys <- xs + eps[,1] > 0
xo1 <- runif(500)
yo1 <- xo1 + eps[,2]
xo2 <- runif(500)
yo2 <- xo2 + eps[,3]
heckit5test <- heckit( ys~xs, list( yo1 ~ xo1, yo2 ~ xo2 ) )
print( heckit5test )
print( summary( heckit5test ) )
print( coef( heckit5test ) )
print( coef( heckit5test, part = "outcome" ) )
print( coef( summary( heckit5test ) ) )
print( coef( summary( heckit5test ), part = "outcome" ) )
print( vcov( heckit5test ) )
print( vcov( heckit5test, part = "outcome" ) )
nobs( heckit5test )
nObs( heckit5test )
fitted( heckit5test, part = "outcome" )
fitted( heckit5test, part = "selection" )
residuals( heckit5test, part = "outcome" )
residuals( heckit5test, part = "selection", type = "response" )
print( model.matrix( heckit5test, part = "outcome" ) )
print( model.matrix( heckit5test, part = "selection" ) )
print( model.frame( heckit5test ) )
