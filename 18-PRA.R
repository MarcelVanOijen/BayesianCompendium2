## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 18. Probabilistic Risk Analysis

  library(mvtnorm)

# Linear and nonlinear datasets for sampling-based single-threshold PRA.
# Dashed lines indicate threshold (black), overall mean (blue), conditional
# expectations (red).
  set.seed(1) ; n <- 3.e2
  mx  <- mz <- 0 ; Vx  <- Vz <- 1 ; rxz <- 0.9
  mxz <- c(mx,mz)
  Sxz <- diag( c(Vx,Vz) ) ; Sxz[1,2] <- Sxz[2,1] <- rxz * sqrt(Vx * Vz)
  xz  <- rmvnorm( n, mxz, Sxz )
  x.l <- xz[,1]  ; z.l <- xz[,2] ; thr.l <- 0 ; H <- which(x.l < thr.l)
  Ez  <- mean( z.l ) ; Ez_H <- mean( z.l[H] ) ; Ez_notH <- mean( z.l[-H] )
  par( mfrow=c(1,2), mar=c(4,4,2,2) )
  plot  ( x.l, z.l, xlab="x", ylab="z", pch=1, col="blue", bty="l", main="Linear dataset" )
  abline( v=thr.l, col="black", lty=2, lwd=2 )
  abline( h=Ez   , col="blue" , lty=2, lwd=2 )
  text  ( x=thr.l   , y=min(z.l), adj=c(-0.2,0.2), font=2, "thr"              )
  text  ( x=max(x.l), y=Ez      , adj=c( 0.8,1.4), font=2, "E[z]", col="blue" )
  lines ( x=c(thr.l,max(x.l)), y=rep(Ez_notH,2), col="red", lty=2, lwd=2 )
  lines ( x=c(min(x.l),thr.l), y=rep(Ez_H,2)   , col="red", lty=2, lwd=2 )
  
  set.seed(1)
  sz     <- 0.1
  x.nl   <- runif( n, 0, 3 ) ; ez <- rnorm( n, 0, sz) ; z.nl <- 1-exp(-2*x.nl) + ez
  thr.nl <- 1 ; H <- which(x.nl < thr.nl)
  Ez     <- mean( z.nl ) ; Ez_H <- mean( z.nl[H] ) ; Ez_notH <- mean( z.nl[-H] )
  plot  ( x.nl, z.nl, xlab="x", ylab="z", col="blue", bty="l", xaxt="n", main="Nonlinear dataset" )
  axis(1, at = 0:3)
  abline( v=thr.nl, col="black", lty=2, lwd=2 )
  abline( h=Ez    , col="blue" , lty=2, lwd=2 )
  text  ( x=thr.nl, y=min(z.nl), adj=c(-0.2,0.2), font=2, "thr"              )
  text  ( x=max(x.nl), y=Ez    , adj=c( 0.8,1.4), font=2, "E[z]", col="blue" )
  lines ( x=c(thr.nl,max(x.nl)), y=rep(Ez_notH,2), col="red", lty=2, lwd=2 )
  lines ( x=c(min(x.nl),thr.nl), y=rep(Ez_H,2)   , col="red", lty=2, lwd=2 )

  PRAn <- function( x, z, thr=-1:1 ) {
    n  <- length(z) ; nthr <- length(thr) ; pH <- V <- R <- rep(NA,nthr)
    Ez_notH <- mean( z[ x>=thr[nthr] ] )
    for(i in 1:nthr) {
      Hi    <- if(i==1) which(x<thr[1]) else which(thr[i-1]<=x & x<thr[i])
      nHi   <- length(Hi) ; Ez_Hi <- mean( z[ Hi ] )
      pH[i] <- nHi / n    ; V[i]  <- Ez_notH - Ez_Hi ; R[i] <- pH[i] * V[i] }
    R.sum <- sum(R) ; pH.sum <- sum(pH) ; V.wsum <- R.sum / pH.sum
    return( list( pH.sum=pH.sum, V.wsum=V.wsum, R.sum=R.sum,
                  pH    =pH    , V     =V     , R    =R     ) ) }
  thrs.l  <- seq(-2  ,0,0.5) ; pran.l  <- PRAn(x.l , z.l , thrs.l )
  thrs.nl <- seq( 0.2,1,0.2) ; pran.nl <- PRAn(x.nl, z.nl, thrs.nl)
  
# Multiple-threshold PRA on two datasets. Five hazardous intervals were
# distinguished; bar labels are interval upper bounds for the environmental
# variable $x$.
  par( mfrow=c(3,2), mar=c(3,3,5,0) )
  barplot( pran.l$pH , xlab="", ylab="" )
    mtext(("p[H]"),side=3,line=1,cex=0.8)
    mtext(~underline("Linear dataset"  ),side=3,line=3)
  barplot( pran.nl$pH, xlab="", ylab="" )
    mtext(("p[H]"),side=3,line=1,cex=0.8)
    mtext(~underline("Nonlinear dataset"),side=3,line=3)
  barplot( pran.l$V  , xlab="", ylab="" )
    mtext(("V"   ),side=3,line=1,cex=0.8)
  barplot( pran.nl$V , xlab="", ylab="" )
    mtext(("V"   ),side=3,line=1,cex=0.8)
  barplot( pran.l$R  , xlab="", ylab="", names.arg=thrs.l  )
    mtext(("R"   ),side=3,line=1,cex=0.8)
  barplot( pran.nl$R , xlab="", ylab="", names.arg=thrs.nl )
    mtext(("R"   ),side=3,line=1,cex=0.8)

  PRAbiGauss <- function( x, z, thr=0 ) {
    mx     <- mean(x) ; sx <- sd(x) ; sz <- sd(z) ; rho <- cor(x,z)
    px.thr <- dnorm(thr,mx,sx)      ; Fx.thr <- pnorm(thr,mx,sx)
    pH     <- Fx.thr
    V      <- rho * sx * sz * px.thr / Fx.thr / (1-Fx.thr)
    R      <- pH * V
    return( c( pH=pH, V=V, R=R ) ) }
  PRAbiGauss( x.l , z.l , thr.l  )
  PRAbiGauss( x.nl, z.nl, thr.nl )

  PRA.UQ <- function( x, z, thr=0 ) {
    n  <- length(z) ; H    <- which(x < thr) ; nH      <- length(H)
    Ez <- mean( z ) ; Ez_H <- mean( z[H] )   ; Ez_notH <- mean( z[-H] )
    pH <- nH / n    ; V    <- Ez_notH - Ez_H ; R       <- Ez_notH - Ez
    a          <- 1 + nH ; b <- 1 + n - nH
    s2_Ez      <- var(z    ) /  n
    s2_Ez_H    <- var(z[ H]) /    nH
    s2_Ez_notH <- var(z[-H]) / (n-nH)  
    s_pH       <- sqrt( a*b / (a+b)^2 / (a+b+1) )
    s_V        <- sqrt( s2_Ez_notH + s2_Ez_H )
    s_R        <- sqrt( s2_Ez_notH + s2_Ez - 2 * s2_Ez_notH * (1-nH/n) )
    return( c( pH=pH, V=V, R=R, s_pH=s_pH, s_V=s_V, s_R=s_R ) ) }
  PRA.UQ( x.l , z.l , thr.l  )
  PRA.UQ( x.nl, z.nl, thr.nl )
