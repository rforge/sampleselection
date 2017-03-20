library( "censReg" )
library( "plm" )

# load outputs that were previously produced by this script 
saved <- new.env()
load( "censRegPanelTest.RData.save", envir = saved )

options( digits = 5 )

printAll <- function( objName ) {
   x <- get( objName )
   if( !exists( objName, envir = saved ) ) {
      cat( "previously saved object '", objName, "' not found\n" )
   } else {
      xSaved <- get( objName, envir = saved )
      if( !isTRUE( all.equal( class( x ), class( xSaved ) ) ) ) {
         cat( "objects '", objName, "' have different classes:\n", sep = "" )
         print( class( x ) )
         print( class( xSaved ) )
      } else if( !isTRUE( all.equal( names( x ), names( xSaved ) ) ) ) {
         cat( "components of objects '", objName, "' have different names:\n",
            sep = "" )
         print( names( x ) )
         print( names( xSaved ) )
      }
      for( n in names( x ) ) {
         if( ! n %in% c( "code", "gradient", "iterations", "last.step",
               "message" ) ) {
            testRes <-  all.equal( x[[ n ]], xSaved[[ n ]], tol = 5e-3 )
            if( !isTRUE( testRes ) ) {
               cat( "component '", n, "' of objects '", objName, "' differ:\n",
                  sep = "" )
               print( testRes )
               print( x[[ n ]] )
               print( xSaved[[ n ]] )
            }
         }
      }
   }
   
   print( x, digits = 1 )
   print( x, logSigma = FALSE , digits = 1 )
   print( maxLik:::summary.maxLik( x ), digits = 1 )
   print( summary( x ), digits = 1 )
   print( summary( x ), logSigma = FALSE , digits = 1 )
   print( round( coef( x ), 2 ) )
   print( round( coef( x, logSigma = FALSE ), 2 ) )
   print( round( vcov( x ), 2 ) )
   print( round( vcov( x, logSigma = FALSE ), 2 ) )
   print( round( coef( summary( x ) ), 2 ) )
   print( round( coef( summary( x ), logSigma = FALSE ), 2 ) )
   try( margEff( x ) )
   print( logLik( x ) )
   print( nobs( x ) )
   print( extractAIC( x ) )
   
   for( n in names( x ) ) {
      cat( "$", n, "\n", sep = "" )
      if( n %in% c( "estimate", "gradientObs" ) ) {
         print( round( x[[ n ]], 2 ) )
      } else if( n %in% c( "hessian" ) ) {
         print( round( x[[ n ]], 1 ) )
      } else if( n %in% c( "gradient" ) ) {
      } else if( ! n %in% c( "last.step" ) ) {
         print( x[[ n ]] )
      }
      cat( "\n" )
   }
   cat( "class\n" )
   print( class( x ) )
}

nId <- 15
nTime <- 4

set.seed( 123 )
pData <- data.frame(
   id = rep( paste( "F", 1:nId, sep = "_" ), each = nTime ),
   time = rep( 1980 + 1:nTime, nId ) )
pData$ui <- rep( rnorm( nId ), each = nTime )
pData$x1 <- rnorm( nId * nTime )
pData$x2 <- runif( nId * nTime )
pData$ys <- -1 + pData$ui + 2 * pData$x1 + 3 * pData$x2 + rnorm( nId * nTime )
pData$y <- ifelse( pData$ys > 0, pData$ys, 0 )
nData <- pData # save data set without information on panel structure
pData <- pdata.frame( pData, c( "id", "time" ) )


## Newton-Raphson method
randEff <- censReg( y ~ x1 + x2, data = pData )
printAll( "randEff" )


## BHHH method
randEffBhhh <- censReg( y ~ x1 + x2, data = pData, method = "BHHH" )
printAll( "randEffBhhh" )


## BFGS method (optim)
randEffBfgs <- censReg( y ~ x1 + x2, data = pData, method = "BFGS" )
printAll( "randEffBfgs" )


## BFGS method (R)
randEffBfgsr <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR" )
printAll( "randEffBfgsr" )


## BHHH with starting values
randEffBhhhStart <- censReg( y ~ x1 + x2, data = pData, method = "BHHH",
   start = c( -0.4, 1.7, 2.2, -0.1, -0.01 ) )
printAll( "randEffBhhhStart" )


## left-censoring at 5
pData$yAdd <- pData$y + 5
randEffAdd <- censReg( yAdd ~ x1 + x2, data = pData, method = "BFGSR", left = 5 )
printAll( "randEffAdd" )


## right-censoring
pData$yNeg <- - pData$y
randEffNeg <- censReg( yNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = 0 )
printAll( "randEffNeg" )


## right-censoring at -5
pData$yAddNeg <- - pData$yAdd
randEffAddNeg <- censReg( yAddNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = -5 )
printAll( "randEffAddNeg" )


## both right and left censoring
pData$yBoth <- ifelse( pData$y < 3, pData$y, 3 )
randEffBoth <- censReg( yBoth ~ x1 + x2, data = pData, method = "BFGSR",
   left = 0, right = 3 )
printAll( "randEffBoth" )


## re-order observations/individuals
set.seed( 234 )
perm <- sample( nId )
nData2 <- nData
nData2$id <- NA
for( i in 1:nId ) {
   nData2$id[ nData$id == paste( "F", i, sep = "_" ) ] <-
      paste( "G", perm[ i ], sep = "_" )
}
pData2 <- pdata.frame( nData2, c( "id", "time" ) )
randEffBfgsr2 <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR" )
all.equal( randEffBfgsr2[ -c(3,5,6,7,9,11,14) ],
   randEffBfgsr[ -c(3,5,6,7,9,11,14) ], tolerance = 1e-2 )

# check if the order of observations/individuals influences the likelihood values
d1c1 <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR", start = coef(randEffBfgsr),
   iterlim = 0 )
all.equal( d1c1[-c(5,6,7,9,12,14,18)], randEffBfgsr[-c(5,6,7,9,12,14,18)] )
d1c1$maximum -  randEffBfgsr$maximum

d2c2 <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR", start = coef(randEffBfgsr2),
   iterlim = 0 )
all.equal( d2c2[-c(5,6,7,9,12,14,18)], randEffBfgsr2[-c(5,6,7,9,12,14,18)] )
d2c2$maximum -  randEffBfgsr2$maximum

d1c2 <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR", 
   start = coef(randEffBfgsr2), iterlim = 0 )
d2c2$maximum - d1c2$maximum
d2c2$gradient - d1c2$gradient

d2c1 <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR", 
   start = coef(randEffBfgsr), iterlim = 0 )
d1c1$maximum - d2c1$maximum
d1c1$gradient - d2c1$gradient

round( d2c2$maximum - d2c1$maximum, 3 )
round( d1c1$maximum - d1c2$maximum, 3 )

d1cS <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR", 
   start = randEffBfgsr$start, iterlim = 0 )
d2cS <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR", 
   start = randEffBfgsr$start, iterlim = 0 )
d1cS$maximum - d2cS$maximum
d1cS$gradient - d2cS$gradient


## unbalanced panel data
nDataUnb <- nData[ -c( 2, 5, 6, 8 ), ]
pDataUnb <- pdata.frame( nDataUnb, c( "id", "time" ) )
randEffBfgsrUnb <- censReg( y ~ x1 + x2, data = pDataUnb, method = "BFGSR" )
printAll( "randEffBfgsrUnb" )


## NAs in data
pDataNa <- pData
obsNa <- which( ! rownames( pData ) %in% rownames( pDataUnb ) )
pDataNa$y[ obsNa[ 1:2 ] ] <- NA
pDataNa$x1[ obsNa[ 3 ] ] <- NA
pDataNa$x2[ obsNa[ c( 1, 2, 4 ) ] ] <- NA
randEffBfgsrNa <- censReg( y ~ x1 + x2, data = pDataNa, method = "BFGSR" )
all.equal( randEffBfgsrNa[ -14 ], randEffBfgsrUnb[ -14 ] )


# returning log-likelihood contributions only (no estimations)
logLikRandEff <- censReg( y ~ x1 + x2, data = pData, start = coef( randEff ),
   logLikOnly = TRUE )
print( logLikRandEff, digits = 1 )
all.equal( sum( logLikRandEff ), c( logLik( randEff ) ) )
logLikStart <- censReg( y ~ x1 + x2, data = pData, 
   start = c( -0.4, 1.7, 2.2, -0.1, -0.01 ), logLikOnly = TRUE )
print( round( c( logLikStart ), 3 ) )
print( round( attr( logLikStart, "gradient" ), 2 ) )


# save all objectives that were produced in this script
# (in order to compare them with objects created in the future)
save.image( "censRegPanelTest.RData" )

