print.aidsEst <- function( x, ... ) {
   cat( "\nDemand analysis with the Almost Ideal " )
   cat( "Demand System (AIDS)\n" )
   cat( "Estimation Method: " )
   if( substr( x$method, 1, 2 ) == "LA" ) {
      cat( "Linear Approximation (LA) with " )
      if( x$px == "S" ) {
         cat( "Stone Index (S)\n" )
      } else if( x$px == "SL" ) {
         cat( "lagged Stone Index (SL)\n" )
      } else if( x$px == "P" ) {
         cat( "Paasche Index (P)\n" )
      } else if( x$px == "L" ) {
         cat( "Laspeyres Index (L)\n" )
      } else if( x$px == "T" ) {
         cat( "Tornqvist Index (T)\n" )
      } else {
         cat( "unknown price index\n" )
      }
   } else if( substr( x$method, 1, 2 ) %in% c( "MK", "IL" ) ) {
      cat( "'Iterated Linear Least Squares Estimator' (IL) starting with " )
      if( x$px == "S" ) {
         cat( "Stone Index (S)\n" )
      } else if( x$px == "SL" ) {
         cat( "lagged Stone Index (SL)\n" )
      } else if( x$px == "P" ) {
         cat( "Paasche Index (P)\n" )
      } else if( x$px == "L" ) {
         cat( "Laspeyres Index (L)\n" )
      } else if( x$px == "T" ) {
         cat( "Tornqvist Index (T)\n" )
      } else {
         cat( "unknown price index\n" )
      }
   }
   cat( "Estimated Coefficients:\n" )
   cat( "alpha\n" )
   print( x$coef$alpha )
   cat( "beta\n" )
   print( x$coef$beta )
   cat( "gamma\n" )
   print( x$coef$gamma )
   if( !is.null( x$coef$delta ) ){
      cat( "delta\n" )
      print( x$coef$delta )
   }
   invisible( x )
}
