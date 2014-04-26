library( "sampleSelection" )
library( "mvtnorm" )
options( digits = 3 )
N <- 1500
NNA <- 5
vc <- diag(3)
vc[lower.tri(vc)] <- c(0.9, 0.5, 0.6)
vc[upper.tri(vc)] <- vc[lower.tri(vc)]
set.seed(1)
## ------- Tobit-5 example ---------
eps <- rmvnorm( N, rep(0, 3), vc )
xs <- runif(N)
ys <- xs + eps[,1] > 0
xo1 <- runif(N)
yo1 <- xo1 + eps[,2]
xo2 <- runif(N)
yo2 <- xo2 + eps[,3]
## Put some NA-s into the data
ys[sample(N, NNA)] <- NA
xs[sample(N, NNA)] <- NA
xo1[sample(N, NNA)] <- NA
xo2[sample(N, NNA)] <- NA
yo1[sample(N, NNA)] <- NA
yo2[sample(N, NNA)] <- NA
testTobit5TwoStep <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), method="2step")
print( testTobit5TwoStep )
print( summary( testTobit5TwoStep ) )
print( coef( testTobit5TwoStep ) )
print( coef( testTobit5TwoStep, part = "outcome" ) )
print( coef( summary( testTobit5TwoStep ) ) )
print( coef( summary( testTobit5TwoStep ), part = "outcome" ) )
print( vcov( testTobit5TwoStep ) )
print( vcov( testTobit5TwoStep, part = "outcome" ) )
nObs( testTobit5TwoStep )
print( fitted( testTobit5TwoStep, part = "outcome" ) )
print( fitted( testTobit5TwoStep, part = "selection" ) )
print( residuals( testTobit5TwoStep, part = "outcome" ) )
print( residuals( testTobit5TwoStep, part = "selection", type = "response" ) )
print( model.matrix( testTobit5TwoStep, part = "outcome" ) )
print( model.matrix( testTobit5TwoStep, part = "selection" ) )
print( model.frame( testTobit5TwoStep ) )
try( logLik( testTobit5TwoStep ) )

testTobit5Ml <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), method="ml")
print( testTobit5Ml )
print( summary( testTobit5Ml ) )
print( coef( testTobit5Ml ) )
print( coef( testTobit5Ml, part = "outcome" ) )
print( coef( summary( testTobit5Ml ) ) )
print( coef( summary( testTobit5Ml ), part = "outcome" ) )
print( vcov( testTobit5Ml ) )
print( vcov( testTobit5Ml, part = "outcome" ) )
nObs( testTobit5Ml )
print( fitted( testTobit5Ml, part = "outcome" ) )
print( fitted( testTobit5Ml, part = "selection" ) )
print( residuals( testTobit5Ml, part = "outcome" ) )
print( residuals( testTobit5Ml, part = "selection" ) )
mmsTestTobit5Ml <- model.matrix( testTobit5Ml, part = "selection" )
print( mmsTestTobit5Ml )
mmoTestTobit5Ml <- model.matrix( testTobit5Ml, part = "outcome" )
print( mmoTestTobit5Ml )
mfTestTobit5Ml <- model.frame( testTobit5Ml )
print( mfTestTobit5Ml )
logLik( testTobit5Ml )

# ML with model.matrices returned
testTobit5MlMm <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ),
   method = "ml", xs = TRUE, xo = TRUE )
mmsTestTobit5MlMm <- model.matrix( testTobit5MlMm, part = "selection" )
attributes( mmsTestTobit5Ml )$assign <- NULL
all.equal( mmsTestTobit5Ml, mmsTestTobit5MlMm )
mmoTestTobit5MlMm <- model.matrix( testTobit5MlMm, part = "outcome" )
attributes( mmoTestTobit5Ml[[ 1 ]] )$assign <- NULL
attributes( mmoTestTobit5Ml[[ 2 ]] )$assign <- NULL
all.equal( mmoTestTobit5Ml, mmoTestTobit5MlMm )
# ML with model.frames returned
testTobit5MlMf <- selection( ys~xs, list( yo1 ~ xo1, yo2 ~ xo2 ),
   method = "ml", mfs = TRUE, mfo = TRUE )
mfTestTobit5MlMf <- model.frame( testTobit5MlMf )
all.equal( mfTestTobit5Ml, mfTestTobit5MlMf )

# return just the model.frame
selection( ys~xs, list( yo1 ~ xo1, yo2 ~ xo2 ), method="model.frame" )

# factors as dependent variable (from Achim Zeileis)
selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ), method = "2step" )
selection( factor( ys ) ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ), method = "2step" )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ), method = "2step" )

selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ) )
selection( factor( ys ) ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ) )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ) )

# with pre-defined list of outcome equations (works since revision 1420)
oList <- list( yo1 ~ xo1, yo2 ~ xo2 )
selection( ys ~ xs, oList, method = "2step" )
selection( factor( ys ) ~ xs, oList, method = "2step" )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs, oList, method = "2step" )

selection( ys ~ xs, oList )
selection( factor( ys ) ~ xs, oList )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs, oList )

# return just the model.frame
selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ), method="model.frame" )
selection( factor( ys ) ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ), method="model.frame" )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ), method="model.frame" )

# the case without intercepts 
cat("Now run tobit5 without intercepts\n")
print(coef(selection( ys ~ xs - 1, list( yo1 ~ xo1 - 1, yo2 ~ xo2 - 1))))
# return just the model.frame
selection( ys ~ xs - 1, list( yo1 ~ xo1 - 1, yo2 ~ xo2 - 1 ),
   method = "model.frame" )

## estimations withs weights that do not work
testTobit5TwoStepWe <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2), 
   method = "2step", weights = rep( 0.5, N ) )
all.equal( testTobit5TwoStepWe[-8], testTobit5TwoStep[-8] )

testTobit5MlWe <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2), 
   method = "ml", weights = rep( 0.5, N ) )
all.equal( testTobit5MlWe[-17], testTobit5Ml[-17] )


## ------- Tobit-2 exmple -----------
vc <- diag(2)
vc[2,1] <- vc[1,2] <- -0.7
eps <- rmvnorm( N, rep(0, 2), vc )
xs <- runif(N)
ys <- xs + eps[,1] > 0
xo <- runif(N)
yo <- (xo + eps[,2])*(ys > 0)
xs[sample(N, NNA)] <- NA
ys[sample(N, NNA)] <- NA
xo[sample(N, NNA)] <- NA
yo[sample(N, NNA)] <- NA
testTobit2TwoStep <- selection(ys~xs, yo ~xo, method="2step")
print( testTobit2TwoStep )
print( summary( testTobit2TwoStep ) )
print( coef( testTobit2TwoStep ) )
print( coef( testTobit2TwoStep, part = "outcome" ) )
print( coef( summary( testTobit2TwoStep ) ) )
print( coef( summary( testTobit2TwoStep ), part = "outcome" ) )
print( vcov( testTobit2TwoStep ) )
print( vcov( testTobit2TwoStep, part = "outcome" ) )
print( testTobit2TwoStep$invMillsRatio )
nObs( testTobit2TwoStep )
print( fitted( testTobit2TwoStep, part = "outcome" ) )
print( fitted( testTobit2TwoStep, part = "selection" ) )
print( residuals( testTobit2TwoStep, part = "outcome" ) )
print( residuals( testTobit2TwoStep, part = "selection", type = "response" ) )
print( model.matrix( testTobit2TwoStep, part = "outcome" ) )
print( model.matrix( testTobit2TwoStep, part = "selection" ) )
print( model.frame( testTobit2TwoStep ) )
try( logLik( testTobit2TwoStep ) )

testTobit2Ml <- selection(ys~xs, yo ~xo, method="ml")
print( testTobit2Ml )
print( summary( testTobit2Ml ) )
print( coef( testTobit2Ml ) )
print( coef( testTobit2Ml, part = "outcome" ) )
print( coef( summary( testTobit2Ml ) ) )
print( coef( summary( testTobit2Ml ), part = "outcome" ) )
print( vcov( testTobit2Ml ) )
print( vcov( testTobit2Ml, part = "outcome" ) )
nObs( testTobit2Ml )
print( fitted( testTobit2Ml, part = "outcome" ) )
print( fitted( testTobit2Ml, part = "selection" ) )
print( residuals( testTobit2Ml, part = "outcome" ) )
print( residuals( testTobit2Ml, part = "selection" ) )
mmsTestTobit2Ml <- model.matrix( testTobit2Ml, part = "selection" )
print( mmsTestTobit2Ml )
mmoTestTobit2Ml <- model.matrix( testTobit2Ml, part = "outcome" )
print( mmsTestTobit2Ml )
mfTestTobit2Ml <- model.frame( testTobit2Ml )
print( mfTestTobit2Ml )
logLik( testTobit2Ml )

# ML with model.matrices returned
testTobit2MlMm <- selection( ys ~ xs, yo ~ xo, method = "ml", 
   xs = TRUE, xo = TRUE )
mmsTestTobit2MlMm <- model.matrix( testTobit2MlMm, part = "selection" )
attributes( mmsTestTobit2Ml )$assign <- NULL
all.equal( mmsTestTobit2Ml, mmsTestTobit2MlMm )
mmoTestTobit2MlMm <- model.matrix( testTobit2MlMm, part = "outcome" )
attributes( mmoTestTobit2Ml )$assign <- NULL
all.equal( mmoTestTobit2Ml, mmoTestTobit2MlMm )
# ML with model.frames returned
testTobit2MlMf <- selection(ys~xs, yo ~xo, method="ml", mfs = TRUE, mfo = TRUE)
mfTestTobit2MlMf <- model.frame( testTobit2MlMf )
all.equal( mfTestTobit2Ml, mfTestTobit2MlMf )
                           # attributes (terms) differ here, I don't exactly know how to improve that

# return just the model.frame
selection( ys~xs, yo ~xo, method = "model.frame" )


## two-step estimation with equal weights
we <- rep( 0.7, N )
testTobit2TwoStepWe <- selection( ys~xs, yo ~xo, method="2step", weights = we )
summary( testTobit2TwoStepWe )
all.equal( coef( testTobit2TwoStepWe ), coef( testTobit2TwoStep ) )
nObs( testTobit2TwoStepWe )
try( logLik( testTobit2TwoStepWe ) )

## ML estimation with equal weights
testTobit2MlWe <- selection( ys~xs, yo ~xo, weights = we )
summary( testTobit2MlWe )
all.equal( coef( testTobit2MlWe ), coef( testTobit2Ml ), tol = 1e-4 )
nObs( testTobit2MlWe )
logLik( testTobit2MlWe )

## two-step estimation with unequal weights
wu <- 2 * runif( N )
testTobit2TwoStepWu <- selection( ys~xs, yo ~xo, method="2step", weights = wu )
summary( testTobit2TwoStepWu )
nObs( testTobit2TwoStepWu )
try( logLik( testTobit2TwoStepWu ) )

## ML estimation with unequal weights
testTobit2MlWu <- selection( ys~xs, yo ~xo, weights = wu )
summary( testTobit2MlWu )
nObs( testTobit2MlWu )
logLik( testTobit2MlWu )

## estimations with weights that do not work
try( selection( ys~xs, yo ~xo, method = "2step", weights = 1:99 ) )

try( selection( ys~xs, yo ~xo, method = "ml", weights = 4:14 ) )


# factors as dependent variable (from Achim Zeileis)
selection( ys ~ xs, yo ~ xo, method = "2step" )
selection( factor( ys ) ~ xs, yo ~ xo, method = "2step" )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs, yo ~ xo,
   method = "2step" )
selection( ys ~ xs, yo ~ xo )
selection( factor( ys ) ~ xs, yo ~ xo )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs, yo ~ xo )
# return just the model.frame
selection( ys ~ xs, yo ~ xo, method = "model.frame" )
selection( factor( ys ) ~ xs, yo ~ xo, method = "model.frame" )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs, yo ~ xo,
   method = "model.frame" )

# the case without intercepts (by Lucas Salazar)
cat("Now run tobit2 without intercepts\n")
print(coef(selection( ys ~ xs - 1, yo ~ xo - 1)))
# return just the model.frame
selection( ys ~ xs - 1, yo ~ xo - 1, method = "model.frame" )

# NA-s in data frames (Nelson Villoria)
set.seed( 98765 )
vc <- diag(2)
vc[2,1] <- vc[1,2] <- -0.8
eps <- rmvnorm( N, rep(0, 2), vc )
xs <- runif(N)
ys <- xs + eps[,1] > 0
xo <- runif(N)
yo <- (xo + eps[,2])*(ys > 0)
xs[sample(N, NNA)] <- NA
ys[sample(N, NNA)] <- NA
xo[sample(N, NNA)] <- NA
yo[sample(N, NNA)] <- NA
data <- data.frame(ys, xs, yo, xo)
rm(eps, xs, ys, xo, yo)
testTobit2ML <- selection(ys~xs, yo ~xo, data=data, method="ml")
print(summary(testTobit2ML))
nObs(testTobit2ML)
logLik(testTobit2ML)

## Raphael Abiry: does 'start' argument work?
init <- coef(testTobit2ML)
testTobit2ML <- selection(ys~xs, yo ~xo, data=data, method="ml", start=init)
print(summary(testTobit2ML))
                           # Note: should be only 1 iteration

## Chris Hane: 'fitted' method and a little complex models
data <- cbind(data, xF=rbinom(nrow(data), 1, 0.5))
testComplex <- selection(ys~xs + factor(xF), yo ~xo, data=data, method="ml")
nObs( testComplex )
f <- fitted(testComplex, "selection")
logLik( testComplex )
