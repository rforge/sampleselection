library( "censReg" )

data( "Affairs", package = "AER" )

estNR <- try( censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[Affairs$affairs > 10,] ) )

estBHHH <- try( censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[ Affairs$affairs > 10, ], 
   method = "BHHH" ) )

estBFGS <- try( censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[ Affairs$affairs > 10, ], 
   method = "BFGS" ) )

estBFGSR <- try( censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[ Affairs$affairs > 10, ], 
   method = "BFGSR" ) )
