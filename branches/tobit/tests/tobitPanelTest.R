library( sampleSelection )
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
pData <- pdata.frame( pData, c( "id", "time" ) )


## Newton-Raphson method
system.time( randEff <- tobit( y ~ x1 + x2, data = pData ) )
summary( randEff )

system.time( randEff2 <- tobit( y ~ x1 + x2, data = pData, vecIndGHQ = FALSE ) )
all.equal( randEff, randEff2 )


## BHHH method
system.time( randEffBhhh <- tobit( y ~ x1 + x2, data = pData,
   method = "BHHH" ) )
summary( randEffBhhh )

system.time( randEffBhhh2 <- tobit( y ~ x1 + x2, data = pData,
   method = "BHHH", vecIndGHQ = FALSE ) )
all.equal( randEffBhhh, randEffBhhh2 )


## BFGS method (optim)
system.time( randEffBfgs <- tobit( y ~ x1 + x2, data = pData,
   method = "BFGS" ) )
summary( randEffBfgs )

system.time( randEffBfgs2 <- tobit( y ~ x1 + x2, data = pData,
   method = "BFGS", vecIndGHQ = FALSE ) )
all.equal( randEffBfgs, randEffBfgs2 )


## BFGS method (R)
system.time( randEffBfgsr <- tobit( y ~ x1 + x2, data = pData,
   method = "BFGSR" ) )
summary( randEffBfgsr )

system.time( randEffBfgsr2 <- tobit( y ~ x1 + x2, data = pData,
   method = "BFGSR", vecIndGHQ = FALSE ) )
all.equal( randEffBfgsr, randEffBfgsr2 )

