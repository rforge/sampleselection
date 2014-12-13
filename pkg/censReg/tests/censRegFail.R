library( "censReg" )

data( "Affairs", package = "AER" )

estNR <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[Affairs$affairs > 10,] )
print( estNR )
try( print( summary( estNR ) ) )


estBHHH <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[ Affairs$affairs > 10, ], 
   method = "BHHH" )
print( estBHHH )
try( print( summary( estBHHH ) ) )


estBFGS <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[ Affairs$affairs > 10, ], 
   method = "BFGS" )
print( estBFGS )
try( print( summary( estBFGS ) ) )

estBFGSR <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[ Affairs$affairs > 10, ], 
   method = "BFGSR" )
print( estBFGSR )
try( print( summary( estBFGSR ) ) )
