library( censReg )

data( "Affairs", package = "AER" )
affairsFormula <- affairs ~ age + yearsmarried + religiousness +
   occupation + rating

## usual tobit estimation
tobitResult <- tobit( affairsFormula, data = Affairs )
print.default( tobitResult )
summary( tobitResult )

## usual tobit estimation, BHHH method
tobitResultBhhh <- tobit( affairsFormula, data = Affairs, method = "BHHH" )
print.default( tobitResultBhhh )
summary( tobitResultBhhh )

## usual tobit estimation, BFGS method
tobitResultBfgs <- tobit( affairsFormula, data = Affairs, method = "BFGS" )
print.default( tobitResultBfgs )
summary( tobitResultBfgs )

## usual tobit estimation, NM method
tobitResultNm <- tobit( affairsFormula, data = Affairs, method = "NM" )
print.default( tobitResultNm )
summary( tobitResultNm )

## usual tobit estimation, SANN method
tobitResultSann <- tobit( affairsFormula, data = Affairs, method = "SANN" )
print.default( tobitResultSann )
summary( tobitResultSann )

## usual tobit estimation with user-defined starting values
tobitResultStart <- tobit( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ) )
print.default( tobitResultStart )
summary( tobitResultStart )

## tobit estimation with left-censoring at 5
Affairs$affairsAdd <- Affairs$affairs + 5
tobitResultAdd <- tobit( affairsAdd ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = 5 )
print.default( tobitResultAdd )
summary( tobitResultAdd )

## tobit estimation with right-censoring
Affairs$affairsNeg <- - Affairs$affairs
tobitResultNeg <- tobit( affairsNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = 0 )
print.default( tobitResultNeg )
summary( tobitResultNeg )

## tobit estimation with right-censoring at -5
Affairs$affairsAddNeg <- - Affairs$affairsAdd
tobitResultAddNeg <- tobit( affairsAddNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = -5 )
print.default( tobitResultAddNeg )
summary( tobitResultAddNeg )

## tobit estimation with left and right censoring
tobitResultBoth <- tobit( affairsFormula, data = Affairs, right = 4 )
print.default( tobitResultBoth )
summary( tobitResultBoth )
