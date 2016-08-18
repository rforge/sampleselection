.mybind<-function(x,y) {
  k<-ncol(x)
  z<-rbind(x,y)
  for(i in 1:k) {
    class(z[,i])<-class(x[,i])
  }
  return( z )
}

.mygetdata <- function(object) {
  fla <- object$call$formula
  cfla <- as.character(fla);
  chfla <- paste(cfla[2],cfla[1],cfla[3],sep="");
  if( length(object$call$data) != 0 ) {
    attach(data.frame(eval(object$call$data)))
  }
  if(grepl(".",chfla,fixed=T)) {
    nm<-names(data.frame(eval(object$call$data)))[-1]
    cnm<-paste(nm,collapse="+")
    cnm<-paste("(",cnm,")",sep="")
    cfla2<-gsub(".",cnm,chfla,fixed=T)
    cfla3<-as.formula(cfla2)
    fla<-cfla3
  }
  av<-all.vars(fla)
  nav<-length(av)
  vav<-data.frame(eval(parse(text=av[1])))
  for (i in 2:nav) {
    vav<-cbind(vav,eval(parse(text=av[i])))
  }
  z<-data.frame(vav)
  names(z)<-av
  detach(data.frame(eval(object$call$data)))
  return( z )
}


.findbeta <- function(object,values) {
  fla<-object$call$formula
  cfla<-as.character(fla)
  chfla<-paste(cfla[2],cfla[1],cfla[3],sep="")
  if (grepl(".",chfla,fixed=T)) {
    nm<-names(data.frame(eval(object$call$data)))[-1]
    cnm<-paste(nm,collapse="+")
    cnm<-paste("(",cnm,")",sep="")
    cfla2<-gsub(".",cnm,chfla,fixed=T)
    fla<-as.formula(cfla2)
  }
  dta <- .mybind(.mygetdata(object),c(0,values))
  lx <- model.matrix(as.formula(fla),data=dta)
  nl <- nrow(lx);
  return( lx[nl,] )
}


.defbeta <- function(object,factors) {
  fla<-object$call$formula
  cfla<-as.character(fla)
  chfla<-paste(cfla[2],cfla[1],cfla[3],sep="")
  if( grepl(".",chfla,fixed=T) ) {
    nm<-names(data.frame(eval(object$call$data)))[-1];
    cnm<-paste(nm,collapse="+");
    cnm<-paste("(",cnm,")",sep="");
    cfla2<-gsub(".",cnm,chfla,fixed=T);
    cfla3<-as.formula(cfla2);
    fla<-cfla3
  }
  dta0 <- .mygetdata(object);
  nd <- ncol(dta0);
  if( missing(factors) ) {
    for( i in 1:nd ) {
      if( class(dta0[,i])%in%c("integer") && length(unique(dta0[,i]))==2) {
        dta0[,i]<-factor(dta0[,i])
      }
    }
  }
  av <- all.vars(fla);
  vls <- rep(0,nd);
  for (i in 2:nd) {
    if (class(dta0[,i])%in%c("factor")) {
      vls[i]<-attr(dta0[,i],"levels")[1]
    } else {
      vls[i]<-mean(dta0[,i])
    }
  }
  dta1<-.mybind(dta0,vls);
  lx2<-model.matrix(as.formula(fla),data=dta1)
  nl<-nrow(lx2);
  lx2[nl,]
}

margEff.censReg <- function( object, calcVCov = TRUE, returnJacobian = FALSE,
  values, ... ) {
  
  ## calculate marginal effects on E[y] at the mean explanatory variables
  allPar <- coef( object, logSigma = FALSE )
  
  # check if the model was estimated with panel data
  isPanel <- "sigmaMu" %in% names( allPar )
  
  ## (not for panel data)
  if( isPanel ) {
    stop( "the margEff() method for objects of class 'censReg'",
          " can not yet be used for panel data models" )
  }
  
  sigma <- allPar[ "sigma" ]
  beta <- allPar[ ! names( allPar ) %in% c( "sigma" ) ]
  if( length( object$xMean ) != length( beta ) ){
    stop( "cannot calculate marginal effects due to an internal error:",
          " please contact the maintainer of this package" )
  }
  
  
  if(missing(values)) {
    x4<-object$xMean 
  } else {
    x4 <- object$xMean
    if((values=="default2")[1]) {
      x3 <- .defbeta(object)
    } else {
      x3<-.findbeta(object,values)
    }
    for(i in 1: length(x4)) {
      x4[i]<-x3[i]
    }
  }
  
  xBeta <- crossprod( x4, beta )
  # class(x3)<-class(object$xMean)
  zRight <- ( object$right - xBeta ) / sigma
  zLeft <- ( object$left - xBeta ) / sigma
  result <- beta[ ! names( beta ) %in% c( "(Intercept)" ) ] * 
    ( pnorm( zRight ) - pnorm( zLeft ) )
  names( result ) <- 
    names( beta )[ ! names( beta ) %in% c( "(Intercept)" ) ]
  
  if( calcVCov || returnJacobian ){
    # compute Jacobian matrix
    jac <- matrix( 0, nrow = length( result ), ncol = length( allPar ) )
    rownames( jac ) <- names( result )
    colnames( jac ) <- names( allPar )
    for( j in names( result ) ) {
      for( k in names( allPar )[ -length( allPar ) ] ) {
        jac[ j, k ] <- 
          ( j == k ) * ( pnorm( zRight ) - pnorm( zLeft ) ) -
          ( beta[ j ] * x4[ k ] / sigma ) *
          ( dnorm( zRight ) - dnorm( zLeft ) );
        # print(beta[ j ] * x4[ k ])
      }
      jac[ j, "sigma"] <- 0
      if( is.finite( object$right ) ) {
        jac[ j, "sigma"] <- jac[ j, "sigma"] - ( beta[ j ] / sigma ) *
          dnorm( zRight ) * zRight
      }
      if( is.finite( object$left ) ) {
        jac[ j, "sigma"] <- jac[ j, "sigma"] + ( beta[ j ] / sigma ) *
          dnorm( zLeft ) * zLeft
      }
    }
    if( calcVCov ) {
      attr( result, "vcov" ) <- 
        jac %*% vcov( object, logSigma = FALSE ) %*% t( jac )
    }
    if( returnJacobian ) {
      attr( result, "jacobian" ) <- jac
    }
  }
  
  # degrees of freedom of the residuals
  attr( result, "df.residual" ) <- object$df.residual
  
  class( result ) <- c( "margEff.censReg", class( result ) )
  
  return( result )
}