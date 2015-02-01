treatReg <- function(selection, outcome,
                      data=sys.frame(sys.parent()),
                      weights = NULL,
                      subset,
                      method="ml",
                      start=NULL,
                      ys=FALSE, xs=FALSE,
                      yo=FALSE, xo=FALSE,
                      mfs=FALSE, mfo=FALSE,
                      print.level=0,
                      ...) {
   ## Heckman-style treatment effect models
   ## selection:   formula
   ##              LHS: must be convertable to two-level factor (e.g. 0-1, 1-2, "A"-"B")
   ##              RHS: ordinary formula as in lm()
   ## outcome:     formula
   ##              should include selection outcome
   ## ys, xs, yo, xo, mfs, mfo: whether to return the response, model matrix or
   ##              the model frame of outcome and selection equation(s)
   ## First the consistency checks
   ## ...          additional arguments for tobit2fit and tobit5fit
   type <- 0
   if(!inherits( selection, "formula" )) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
   if(inherits(outcome, "formula")) {
      if( length( outcome ) != 3 ) {
         stop( "argument 'outcome' must be a 2-sided formula" )
      }
   }
   else
       stop("argument 'outcome' must be a formula" )
   if(!missing(data)) {
      if(!inherits(data, "environment") & !inherits(data, "data.frame") & !inherits(data, "list")) {
         stop("'data' must be either environment, data.frame, or list (currently a ", class(data), ")")
      }
   }
   if(print.level > 0)
       cat("Treatment effect model", type, "model\n")
   probitEndogenous <- model.frame( selection, data = data)[ , 1 ]
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   if( !is.null( weights ) && type != 2 ) {
      warning( "argument 'weights' is ignored" )
      weights <- NULL
   }
   ## now check whether two-step method was requested
   cl <- match.call()
   if(method == "2step") {
      twoStep <- heckitTfit(selection, outcome, data=data,
                            weights = weights,
                            print.level = print.level, ... )
      twoStep$call <- cl
      class(twoStep) <- c("selection", class(twoStep))
      return(twoStep)
   }
   ## Now extract model frames etc
   ## YS (selection equation)
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$drop.unused.levels <- TRUE
   mfS$na.action <- na.pass
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mfS <- eval(mfS, parent.frame())
   mtS <- attr(mfS, "terms")
   XS <- model.matrix(mtS, mfS)
   YS <- model.response(mfS)
   YSLevels <- levels( as.factor( YS ) )
   if( length( YSLevels ) != 2 ) {
      stop( "the left hand side of the 'selection' formula has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   YS <- as.integer(YS == YSLevels[ 2 ])
                                        # selection will be kept as integer internally
   ## check for NA-s.  Because we have to find NA-s in several frames, we cannot use the standard na.
   ## functions here.  Find bad rows and remove them later.
   ## We check XS and YS separately, because mfS may be a data frame with complex structure (e.g.
   ## including matrices)
   badRow <- is.na(YS)
   badRow <- badRow | apply(XS, 1, function(v) any(is.na(v)))
   ## YO (outcome equation)
   ## Here we should include a possibility for the user to
   ## specify the model.  Currently just a guess.
   binaryOutcome <- FALSE
   oArg <- match("outcome", names(mf), 0)
                           # find the outcome argument
   m <- match(c("outcome", "data", "subset",
                "offset"), names(mf), 0)
   ## replace the outcome list by the first equation and evaluate it
   mfO <- mf[c(1, m)]
   mfO$drop.unused.levels <- TRUE
   mfO$na.action <- na.pass
   mfO[[1]] <- as.name("model.frame")
                           # eval it as model frame
   names(mfO)[2] <- "formula"
   mfO <- eval(mfO, parent.frame())
                           # Note: if unobserved variables are
                           # marked as NA, eval returns a
                           # subframe of visible variables only.
                           # We have to check it later
   mtO <- attr(mfO, "terms")
   XO <- model.matrix(mtO, mfO)
   YO <- model.response(mfO)
   if(is.logical(YO) |
      (is.factor(YO) & length(levels(YO)) == 2)) {
      binaryOutcome <- TRUE
   }
   badRow <- badRow | is.na(YO)
   badRow <- badRow | (apply(XO, 1, function(v)
      any(is.na(v))))
                           # rows in outcome, which contain 
   if( !is.null( weights ) ) {
      if( length( weights ) != length( badRow ) ) {
         stop( "number of weights (", length( weights ), ") is not equal",
              " to the number of observations (", length( badRow ), ")" )
      }
      badRow <- badRow | is.na( weights )
   }   
   if(print.level > 0) {
      cat(sum(badRow), "invalid observations\n")
   }
   if( method == "model.frame" ) {
      mf <- mfS
      mf <- cbind( mf, mfO[ , ! names( mfO ) %in% names( mf ), drop = FALSE ] )
      return( mf[ !badRow, ] )
   }
   XS <- XS[!badRow,, drop=FALSE]
   YS <- YS[!badRow]
   XO <- XO[!badRow,, drop=FALSE]
   YO <- YO[!badRow]
   weightsNoNA <- weights[ !badRow ]
   YO[ YS == 0 ] <- NA
   XO[ YS == 0, ] <- NA
   NXS <- ncol(XS)
   NXO <- ncol(XO)
   iBetaS <- 1:NXS
   iBetaO <- max(iBetaS) + seq(length=NXO)
   if(!binaryOutcome) {
      iSigma <- max(iBetaO) + 1
      iRho <- max(iSigma) + 1
   }
   else
      iRho <- max(iBetaO) + 1
   nParam <- iRho
   twoStep <- NULL
   if(is.null(start)) {
                           # start values by Heckman 2-step method
      start <- numeric(nParam)
      twoStep <- heckit2fit(selection, outcome, data=data,
                            print.level = print.level, weights = weights )
      coefs <- coef(twoStep, part="full")
      start[iBetaS] <- coefs[twoStep$param$index$betaS]
      if(!binaryOutcome) {
         start[iBetaO] <- coefs[twoStep$param$index$betaO]
         start[iSigma] <- coefs[twoStep$param$index$sigma]
      }
      else
         start[iBetaO] <- coefs[twoStep$param$index$betaO]/coefs[twoStep$param$index$sigma]
      start[iRho] <- coefs[twoStep$param$index$rho]
      if(start[iRho] > 0.99)
         start[iRho] <- 0.99
      else if(start[iRho] < -0.99)
         start[iRho] <- -0.99
   }
   if(is.null(names(start))) {
      if(!binaryOutcome) {
         names(start) <- c(colnames(XS), colnames(XO), "sigma",
                           "rho")
      }
      else
         names(start) <- c(colnames(XS), colnames(XO), 
                           "rho")
   }                                        # add names to start values if not present
   if(!binaryOutcome) {
      estimation <- tobitTfit(YS, XS, YO, XO, start,
                              weights = weightsNoNA,
                              print.level=print.level, ...)
      iErrTerms <- c(sigma=iSigma, rho=iRho )
   }
   else {
      ## estimation <- tobitTBfit(YS, XS, YO, XO, start, weights = weightsNoNA,
      ##                          print.level=print.level, ...)
      ## iErrTerms <- c(rho=iRho)
      stop("Binary outcome models are not implemented")
   }
   param <- list(index=list(betaS=iBetaS,
                 betaO=iBetaO,
                 errTerms=iErrTerms,
                 outcome = iBetaO),
                 NXS=ncol(XS), NXO=ncol(XO),
                 N0=sum(YS==0), N1=sum(YS==1),
                 nObs=length(YS), nParam=length(start),
                 df=length(YS) - length(start),
                 levels=YSLevels
                           # levels[1]: selection 1; levels[2]: selection 2
                 )
   result <- c(estimation,
               twoStep=list(twoStep),
               start=list(start),
               param=list(param),
               call=cl,
               termsS=mtS,
               termsO=mtO,
               ys=switch(as.character(ys), "TRUE"=list(YS), "FALSE"=NULL),
               xs=switch(as.character(xs), "TRUE"=list(XS), "FALSE"=NULL),
               yo=switch(as.character(yo), "TRUE"=list(YO), "FALSE"=NULL),
               xo=switch(as.character(xo), "TRUE"=list(XO), "FALSE"=NULL),
               mfs=switch(as.character(mfs), "TRUE"=list(mfS[!badRow,]), "FALSE"=NULL),
               mfo=switch(as.character(mfs),
               "TRUE"=list(mfO[!badRow,]), "FALSE"=NULL)
               )
   result$binaryOutcome <- binaryOutcome
   class( result ) <- class( estimation ) 
   return(result)
}
