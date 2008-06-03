path <- "sampleSelection/R/"
sampleSelectionFiles <- list.files( path, all.files = TRUE )
sampleSelectionFiles <- grep( ".*\\.R$", sampleSelectionFiles, value = TRUE )
for( i in seq( along = sampleSelectionFiles ) ) {
   returnedByTry <- try( source(
      paste( path, sampleSelectionFiles[ i ], sep = "" ) ) )
   if( class( returnedByTry ) == "try-error" ) {
      cat( paste( sampleSelectionFiles[ i ], ": ", returnedByTry, sep = "" ) )
   }
}
