
options(digits=5)

## Leeman Lucas (and many others): binary outcome

set.seed(0)
N <- 500
rho <- 0.7
library( "mvtnorm" )
eps <- rmvnorm(N, c(0,0), matrix(c(1,rho,rho,1), 2, 2) )
xs <- runif(N)
ysX <- 3*xs + eps[,1]
ys <- ysX > 0
xo <- runif(N)
yoX <- -1 + 2*xo + eps[,2]
yo <- factor((yoX > 0)*(ys > 0))
                           # binary outcome, only observable if ys>0
print(table(ys, yo, exclude=NULL))
library( "sampleSelection" )

# estimation with BHHH method
ss <- selection( ys ~ xs, yo ~ xo, steptol = 1e-12 )
print(summary(ss))

# estimation with BFGS method
ssBFGS <- selection( ys ~ xs, yo ~ xo, maxMethod = "BFGS" )
print( summary( ssBFGS ) )
all.equal( coef( ssBFGS ), coef( ss ), tol = 1e-2 )

# BHHH estimation with equal weights
we <- rep( 0.7, N )
ssWe <- selection( ys ~ xs, yo ~ xo, weights = we, steptol = 1e-12 )
summary( ssWe )
all.equal( coef( ssWe ), coef( ss ), tol = 1e-2 )

# BHHH estimation with equal weights
ssWeBFGS <- selection( ys ~ xs, yo ~ xo, weights = we, maxMethod = "BFGS" )
summary( ssWeBFGS )
all.equal( coef( ssWeBFGS ), coef( ssBFGS ), tol = 1e-2 )

# BHHH estimation with unequal weights
wu <- 2 * runif( N )
ssWu <- selection( ys ~ xs, yo ~ xo, weights = wu )
summary( ssWu )

# BFGS estimation with unequal weights
ssWuBFGS <- selection( ys ~ xs, yo ~ xo, weights = wu, maxMethod = "BFGS" )
summary( ssWuBFGS )
all.equal( coef( ssWuBFGS ), coef( ssWu ), tol = 1e-2 )

# comparison of estimated coefficients, standard errors, and logLik values
round( rbind( coef( ss ), coef( ssBFGS ), coef( ssWe ), coef( ssWeBFGS ),
   coef( ssWu ), coef( ssWuBFGS ) ), 3 )
round( rbind( coef( summary( ss ) )[ , 2 ], coef( summary( ssBFGS ) )[ , 2 ],
   coef( summary( ssWe ) )[ , 2 ], coef( summary( ssWeBFGS ) )[ , 2 ],
   coef( summary( ssWu ) )[ , 2 ], coef( summary( ssWuBFGS ) )[ , 2 ] ), 3 )
print( rbind( logLik( ss ), logLik( ssBFGS ), logLik( ssWe ),
   logLik( ssWeBFGS ), logLik( ssWu ), logLik( ssWuBFGS ) ), digits = 6 )

# binary outcome NA if unobserved
yo[ !ys ] <- NA
print(table(ys, yo, exclude=NULL))
ssN <- selection( ys ~ xs, yo ~ xo, steptol = 1e-12 )
print(summary(ssN))
all.equal(ss,ssN)

# binary outcome logical
yo <- yoX > 0 & ys
print(table(ys, yo, exclude=NULL))
ssL <- selection( ys ~ xs, yo ~ xo, steptol = 1e-12 )
print(summary(ssL))
all.equal(ss,ssL)

# binary outcome logical and NA if unobserved
yo[ !ys ] <- NA
print(table(ys, yo, exclude=NULL))
ssLN <- selection( ys ~ xs, yo ~ xo, steptol = 1e-12 )
print(summary(ssLN))
all.equal(ssL,ssLN)
