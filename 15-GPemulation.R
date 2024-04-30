## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 15. Gaussian Processes and Model Emulation

  library(geoR)

  x <- c(10,20,30) ; y <- c(6.09,8.81,10.66)

  Vy <- 3 ; phi <- 1e1       # GP variance and correlation length
  dx <- as.matrix( dist(x) ) # Distance matrix
  Sy <- exp( -dx/phi ) * Vy  # Covariance matrix

  mb <- c( 0, 0 ) ; Sb <- diag( 1e4, 2 )

  GP.est <- function( x, y, Sy, mb, Sb, X=cbind(1,x) ) {
    Sb_y <- solve( solve(Sb) + t(X) %*% solve(Sy) %*% X )
    mb_y <- Sb_y %*% ( solve(Sb) %*% mb + t(X) %*% solve(Sy) %*% y )
    return( list( "mb_y"=mb_y, "Sb_y"=Sb_y ) ) }

  b_y <- GP.est(x,y,Sy,mb,Sb) ; mb_y <- b_y$mb_y; Sb_y <- b_y$Sb_y
# ( 3.906838 , 3.906838 )
# [ 6.7424941, -0.25922418
#  -0.2592242,  0.01296323 ]

  GP.pred <- function(x0,x,y,Sy,phi,mb_y,Sb_y,X0=c(1,x0),X=cbind(1,x)) {
    dx0  <- if( is.vector(x) ) { abs(x-x0)
            } else { sapply(1:length(y),function(i){dist(rbind(x0,x[i,]))}) }
    C0   <- Sy[1] * exp( -dx0/phi )
    m0_y <- X0 %*% mb_y - t(C0) %*% solve(Sy) %*% (X %*% mb_y - y)
       a <- X0 - t(C0) %*% solve(Sy) %*% X
    S0_y <- Sy[1] - t(C0) %*% solve(Sy) %*% C0 + a %*% Sb_y %*% t(a)
    return( list( "m0_y"=m0_y, "S0_y"=S0_y ) ) }

  x0a <- 10.1 ; x0b <- 15
  g0a <- GP.pred(x0a,x,y,Sy,phi,mb_y,Sb_y) ; m0a_y <- g0a$m0_y ; S0a_y <- g0a$S0_y
  g0b <- GP.pred(x0b,x,y,Sy,phi,mb_y,Sb_y) ; m0b_y <- g0b$m0_y ; S0b_y <- g0b$S0_y
  cat("x0=10.1: mean =",m0a_y,"; variance =",S0a_y) # 6.117023, 0.05926102
  cat("x0=15.0: mean =",m0b_y,"; variance =",S0b_y) # 7.437081, 1.410466

  s.11     <- as.geodata( cbind( x, rep(0,3), y) )
  xpred    <- c(10.1,15) ; xpred.11 <- expand.grid( xpred, 0 )
  model.11 <- model.control(cov.m  ="exponential",
                            trend.d=~coords[,1], trend.l=~ xpred.11[,1])
  prior.11 <- prior.control(beta.prior     ="normal", beta     =mb,
                            beta.var.std=Sb/Vy,
                            sigmasq.prior  ="fixed" , sigmasq  =Vy,
                            phi.prior      ="fixed" , phi      =phi,
                            tausq.rel.prior="fixed" , tausq.rel=0)
  out.11   <- output.control(messages=FALSE)
  geoR.11  <- krige.bayes(s.11,loc=xpred.11,mod=model.11,pr=prior.11,out=out.11)

  with(geoR.11$posterior$beta$pars,{cat("Mean:",mean,"\nVariance:\n");print(var)})
# For phi=10: mean=(3.906838,0.228601), vars=(6.7424941,0.01296323),cov=-0.2592242

  with(geoR.11$predictive,cat("Means:",mean,", Variances:",variance))
# For phi=10: means = 6.117023, 7.437081, Variances = 0.05926102, 1.410466

# GP-emulation of a 1D model y=f(x) by means of geoR: results for 11 different
# covariance functions."}
  xpred <- seq(5,35,l=121) ; xpred.11 <- expand.grid( xpred, 0 )
  covfunctions <- c("exponential", "gaussian", "matern", "spherical", "circular",
                    "cubic", "wave", "powered.exponential", "cauchy",
                    "gneiting", "pure.nugget")
  ncf <- length(covfunctions) ; nr <- round(sqrt(ncf)) ; nc <- ceiling((ncf)/nr)
  par( mfrow=c(nr,nc), omi=c(0,0,0.5,0), mar=c(2, 2, 2, 1) )
  for( cf in covfunctions ) {
    model.11 <- model.control( cov.m  =cf, trend.d=~coords[,1], trend.l=~ xpred.11[,1] )
    geoR.11 <- krige.bayes( s.11, loc=xpred.11, model=model.11,
                            prior=prior.11, output=out.11 )
    ypred   <- geoR.11$predictive$mean ; sdpred <- sqrt(geoR.11$predictive$variance)
    plotrange <- range( 0, ypred + 2*sdpred )
    plot( xpred, ypred, type="l", ylim=plotrange, xlab="x", ylab="E[y] ± 2 s.d.",
          main=cf )
    points( xpred, ypred-2*sdpred, type="l", lty="dashed" )
    points( xpred, ypred+2*sdpred, type="l", lty="dashed" )
    points( x, y, pch=16, col="blue" )
  }

  EXPOL6s <- function( t=0, b=c(1,1,1,1,1,1) ) {
    t      <- 100     * t
    I0     <-  10     * b[1] # MJ PAR m-2 ground d-1
    K      <-   1     * b[2] # m2 ground m-2 leaf area
    LAIMAX <-   3     * b[3] # m2 leaf area m-2 ground
    LAR    <-   0.007 * b[4] # m2 leaf area g-1 DM
    LUE    <-   2     * b[5] # g DM MJ-1 PAR
    W0     <-   1     * b[6] # g DM m-2 ground
    e1     <- exp(I0*K*LAR*LUE*t) ; e2 = exp(K*LAR*W0) ; e3 = exp(-K*LAIMAX)
    W      <- W0*e3 + ( log(e2+1/e1-1)/(K*LAR) + I0*LUE*t ) * (1-e3)
    LAI    <- log( (1+e1*(e2-1)) / (1+e1*(e2-1)*e3) ) / K
    return( list( "W"=W, "LAI"=LAI ) ) }

  bmin6  <- c( 0.99, 0.5, 1/3, 5/7, 1/2, 0.5 ) ; nb <- length(bmin6)
  bmax6  <- c( 1.01, 1.5, 5/3, 9/7, 3/2, 1.5 )
  mb6    <- rowMeans( cbind(bmin6,bmax6) )
  range6 <- bmax6 - bmin6
  var6   <- (range6^2) / 12
  Sb6    <- diag( var6 )

  set.seed(13)
  nl        <- 9 ; lhs6 <- lhs::randomLHS( n=nl, k=nb )
  sample.b  <- t(sapply( 1:nl, function(i){bmin6 + lhs6[i,]*(bmax6-bmin6)} ))
  sample.x  <- c(20,40,60,80,100) / 100 ; nx <- length(sample.x)
  sample.xb <- NULL
  for (x in sample.x) {sample.xb <- rbind( sample.xb, cbind(x,sample.b) ) }

# Training set for GP-emulation of EXPOL6: Outputs from 45 input vectors.
  sample.out <- t(sapply(1:(nl*nx), function(i){unlist(
    EXPOL6s(t=sample.xb[i,1],b=sample.xb[i,2:7]) )} ))
  # head( sample.out[,"W"] ) ; head( sample.out[,"LAI"] )
  par( mfrow=c(1,1), mar=c(4,4,2,2) )
  plot(sample.out,xlab="W (g m-2)",ylab="LAI (m2 m-2)",main="Training set (EXPOL6 outputs)")

  x    <- sample.xb ; dx <- as.matrix( dist(x) )
  y    <- sample.out[,"LAI"]
  Vy   <- 3 ; phi <- 0.2
  Sy   <- exp( -dx/phi ) * Vy
  X    <- cbind(1,x) ; mb <- rep(0,8) ; Sb <- diag(1e4,8)
  b_y  <- GP.est( x, y, Sy, mb=mb, Sb=Sb, X=X )
  mb_y <- b_y$mb_y; Sb_y <- b_y$Sb_y
  cat( "mb_y = (", signif(t(mb_y),4), ")" )
  cat( "Standard deviations: ", signif(sqrt(diag(Sb_y)),4) )

# In-sample emulation of EXPOL6: interpolation of LAI (m$^2$ m$^{-2}$) vs.
# scaled time for the 9 parameter vectors that were used in GP calibration.
# Blue dots are the training points, i.e. the LAI-output values from EXPOL6.
  nr <- ceiling(sqrt(nl)) ; nc <- ceiling((nl)/nr)
  par( mfrow=c(nr,nc), omi=c(0,0,0.5,0), mar=c(2, 2, 2, 1) )
  for( il in 1:nl ) {
    b.pred  <- sample.b[il,] ; nt <- 101 ; t.pred  <- seq(0.1,1.1,l=nt)
    tb.pred <- matrix( NA, nrow=nt, ncol=1+nb )
    for (it in 1:nt) tb.pred[it,] <- c(t.pred[it],b.pred)
    out.pred <- matrix( NA, nrow=nt, ncol=3 )
    for (it in 1:nt) {
      x0    <- tb.pred[it,] ; X0 <- c(1,x0)
      out.p <- unlist( GP.pred( x0, x, y, Sy, phi, mb_y, Sb_y, X0, X ) )
      out.pred[it,] <- c( t.pred[it], out.p[1], sqrt(max(0,out.p[2])) )
    }
    xpred <- out.pred[,1] ; ypred <- out.pred[,2] ; sdpred <- out.pred[,3]
    ib <- which( apply( sample.xb[,2:7], 1, function(b){all(b == b.pred)} ) )
    x.training  <- sample.xb[ib,1]
    y.training  <- sample.out[ib,][,"LAI"]
  
    plotrange <- range( 0, ypred + 2*sdpred )
    plot( xpred, ypred, type="l", ylim=plotrange,
          main=paste0("b[",il,"]"), xlab="t", ylab="E[y] ± 2 s.d." )
    points( xpred, ypred-2*sdpred, type="l", lty="dashed" )
    points( xpred, ypred+2*sdpred, type="l", lty="dashed" )
    points( x.training, y.training, pch=16, col="blue" )
  }

# Out-of-sample emulation of EXPOL6: interpolation of LAI (m$^2$ m$^{-2}$) vs.
# scaled time for 6 parameter vectors that were not used in GP calibration.
# Blue line: EXPOL6. Black line and dashed lines: emulator.
  ntest <- 6
  sample.b  <- t(sapply( 1:ntest, function(i){bmin6 + runif(6)*(bmax6-bmin6)} ))
  par( mfrow=c(2,3), omi=c(0,0,0.5,0), mar=c(2, 2, 2, 1) )
  for( ip in 1:ntest ) {
    b.pred  <- sample.b[ip,]
    nt      <- 101 ; t.pred  <- seq(0.1,1.1,l=nt)
    tb.pred <- matrix( NA, nrow=nt, ncol=1+nb )
    for (it in 1:nt) tb.pred[it,] <- c(t.pred[it],b.pred)
    out.pred <- matrix( NA, nrow=nt, ncol=3 )
    for (it in 1:nt) {
      x0    <- tb.pred[it,] ; X0 <- c(1,x0)
      out.p <- unlist( GP.pred( x0, x, y, Sy, phi, mb_y, Sb_y, X0, X ) )
      out.pred[it,] <- c( t.pred[it], out.p[1], sqrt(max(0,out.p[2])) )
    }
    xpred <- out.pred[,1] ; ypred <- out.pred[,2] ; sdpred <- out.pred[,3]
    plotrange <- range( -4, 8 )
    plot( xpred, ypred, type="l", ylim=plotrange,
          main=paste0("b[",ip,"]"), xlab="t", ylab="E[y] ± 2 s.d." )
    points( xpred, ypred-2*sdpred, type="l", lty="dashed" )
    points( xpred, ypred+2*sdpred, type="l", lty="dashed" )
    
    points( xpred, EXPOL6s(xpred,b.pred)$LAI, type="l", col="blue" )
  }

# Out-of-sample emulation of EXPOL6: estimation of LAI (m$^2$ m$^{-2}$) at 4
# scaled times for 1000 parameter vectors that were not used in GP calibration.
# x-axis: emulator, y-axis: EXPOL6.
  ntest <- 1e3 ; lhs6 <- lhs::randomLHS( n=ntest, k=nb )
  sample.b  <- t(sapply( 1:ntest, function(i){bmin6 + lhs6[i,]*(bmax6-bmin6)} ))
  par( mfrow=c(2,2), mar=c(2,2,2,2) )
  for (ttest in c(0.1,0.4,0.7,1.0)) { 
    sample.pred <- matrix( NA, nrow=ntest, ncol=nb+2 )
    for( ip in 1:ntest ) {
      b.pred  <- sample.b[ip,]
      x0      <- c( ttest, b.pred ) ; X0 <- c(1,x0)
      out.em  <- unlist( GP.pred( x0, x, y, Sy, phi, mb_y, Sb_y, X0, X ) )
      out.mod <- EXPOL6s( ttest, b.pred )$LAI
      sample.pred[ip,] <- c( b.pred, out.em[1], out.mod ) }
    plot( sample.pred[,7:8], main=paste0("t = ",ttest), pch=".",
          xlab="LAI (emulator)", ylab="LAI (model)" )
    sample.cor <- cor( sample.pred )[1:6,7:8]
  }
