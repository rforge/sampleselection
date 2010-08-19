library( censReg )

data( "Affairs", package = "AER" )
affairsFormula <- affairs ~ age + yearsmarried + religiousness +
   occupation + rating

## usual tobit estimation
estResult <- censReg( affairsFormula, data = Affairs )
print.default( estResult )
summary( estResult )
coef( estResult )
coef( estResult, logSigma = FALSE )

## usual tobit estimation, BHHH method
estResultBhhh <- censReg( affairsFormula, data = Affairs, method = "BHHH" )
print.default( estResultBhhh )
summary( estResultBhhh )

## usual tobit estimation, BFGS method
estResultBfgs <- censReg( affairsFormula, data = Affairs, method = "BFGS" )
print.default( estResultBfgs )
summary( estResultBfgs )

## usual tobit estimation, NM method
estResultNm <- censReg( affairsFormula, data = Affairs, method = "NM" )
print.default( estResultNm )
summary( estResultNm )

## usual tobit estimation, SANN method
estResultSann <- censReg( affairsFormula, data = Affairs, method = "SANN" )
print.default( estResultSann )
summary( estResultSann )

## usual tobit estimation with user-defined starting values
estResultStart <- censReg( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ) )
print.default( estResultStart )
summary( estResultStart )

## estimation with left-censoring at 5
Affairs$affairsAdd <- Affairs$affairs + 5
estResultAdd <- censReg( affairsAdd ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = 5 )
print.default( estResultAdd )
summary( estResultAdd )
coef( estResultAdd )
coef( estResultAdd, logSigma = FALSE )

## estimation with right-censoring
Affairs$affairsNeg <- - Affairs$affairs
estResultNeg <- censReg( affairsNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = 0 )
print.default( estResultNeg )
summary( estResultNeg )
coef( estResultNeg )
coef( estResultNeg, logSigma = FALSE )

## estimation with right-censoring at -5
Affairs$affairsAddNeg <- - Affairs$affairsAdd
estResultAddNeg <- censReg( affairsAddNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = -5 )
print.default( estResultAddNeg )
summary( estResultAddNeg )
coef( estResultAddNeg )
coef( estResultAddNeg, logSigma = FALSE )

## estimation with left and right censoring
estResultBoth <- censReg( affairsFormula, data = Affairs, right = 4 )
print.default( estResultBoth )
summary( estResultBoth )
coef( estResultBoth )
coef( estResultBoth, logSigma = FALSE )
