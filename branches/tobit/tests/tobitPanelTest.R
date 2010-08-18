library( censReg )
library( plm )

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
system.time( randEff <- tobit( y ~ x1 + x2, data = pData ) )
summary( randEff )
print.default( randEff )


## BHHH method
system.time( randEffBhhh <- tobit( y ~ x1 + x2, data = pData,
   method = "BHHH" ) )
summary( randEffBhhh )
print.default( randEffBhhh )


## BFGS method (optim)
system.time( randEffBfgs <- tobit( y ~ x1 + x2, data = pData,
   method = "BFGS" ) )
summary( randEffBfgs )
print.default( randEffBfgs )


## BFGS method (R)
system.time( randEffBfgsr <- tobit( y ~ x1 + x2, data = pData,
   method = "BFGSR" ) )
summary( randEffBfgsr )
print.default( randEffBfgsr )


## left-censoring at 5
pData$yAdd <- pData$y + 5
randEffAdd <- tobit( yAdd ~ x1 + x2, data = pData, method = "BFGSR", left = 5 )
summary( randEffAdd )
print.default( randEffAdd )


## right-censoring
pData$yNeg <- - pData$y
randEffNeg <- tobit( yNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = 0 )
summary( randEffNeg )
print.default( randEffNeg )


## right-censoring at -5
pData$yAddNeg <- - pData$yAdd
randEffAddNeg <- tobit( yAddNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = -5 )
summary( randEffAddNeg )
print.default( randEffAddNeg )


## both right and left censoring
pData$yBoth <- ifelse( pData$y < 3, pData$y, 3 )
randEffBoth <- tobit( yBoth ~ x1 + x2, data = pData, method = "BFGSR",
   left = 0, right = 3 )
summary( randEffBoth )
print.default( randEffBoth )


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
system.time( randEffBfgsr2 <- tobit( y ~ x1 + x2, data = pData2,
   method = "BFGSR" ) )
all.equal( randEffBfgsr2[ -11 ], randEffBfgsr[ -11 ] )
all.equal( sort( randEffBfgsr2[[ 11 ]] ), sort( randEffBfgsr[[ 11 ]] ) )


## unbalanced panel data
nDataUnb <- nData[ -c( 2, 5, 6, 8 ), ]
pDataUnb <- pdata.frame( nDataUnb, c( "id", "time" ) )
system.time( randEffBfgsrUnb <- tobit( y ~ x1 + x2, data = pDataUnb,
   method = "BFGSR" ) )
summary( randEffBfgsrUnb )
print.default( randEffBfgsrUnb )


## NAs in data
pDataNa <- pData
obsNa <- which( ! rownames( pData ) %in% rownames( pDataUnb ) )
pDataNa$y[ obsNa[ 1:2 ] ] <- NA
pDataNa$x1[ obsNa[ 3 ] ] <- NA
pDataNa$x2[ obsNa[ c( 1, 2, 4 ) ] ] <- NA
system.time( randEffBfgsrNa <- tobit( y ~ x1 + x2, data = pDataNa,
   method = "BFGSR" ) )
all.equal( randEffBfgsrNa, randEffBfgsrUnb )


