## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 23. Spatial Modelling and Scaling Error

  library(geoR)
  library(mvtnorm)
  library(rsvg)

  GP.est <- function( x, y, Sy, mb, Sb, X=cbind(1,x) ) {
    Sb_y <- solve( solve(Sb) + t(X) %*% solve(Sy) %*% X )
    mb_y <- Sb_y %*% ( solve(Sb) %*% mb + t(X) %*% solve(Sy) %*% y )
    return( list( "mb_y"=mb_y, "Sb_y"=Sb_y ) )
  }

  GP.pred <- function(x0,x,y,Sy,phi,mb_y,Sb_y,X0=c(1,x0),X=cbind(1,x)) {
    dx0  <- if( is.vector(x) ) { abs(x-x0)
            } else { sapply(1:length(y),function(i){dist(rbind(x0,x[i,]))}) }
    C0   <- Sy[1] * exp( -dx0/phi )
    m0_y <- X0 %*% mb_y - t(C0) %*% solve(Sy) %*% (X %*% mb_y - y)
       a <- X0 - t(C0) %*% solve(Sy) %*% X
    S0_y <- Sy[1] - t(C0) %*% solve(Sy) %*% C0 + a %*% Sb_y %*% t(a)
    return( list( "m0_y"=m0_y, "S0_y"=S0_y ) )
  }

  VcondR <- function( W ) {
    n  <- dim(W)[1] ; Vcond <- rep( NA, n )
    Wk <- W         ; R     <- matrix( 0, nrow=n, ncol=n )
    for(k in n:2){
      ik       <- 1 : (k-1)
      Vcond[k] <- 1 / Wk[ k, k]
      R[k,ik]  <-    -Wk[ k,ik] * Vcond[k]
      Wk       <-     Wk[ik,ik] - as.matrix(R[k,ik]) %*% R[k,ik] / Vcond[k]
    }
    Vcond[1] <- 1 / Wk[1,1]
    return( list( Vcond=Vcond, R=zapsmall(R) ) )
  }

  GaussCond <- function( mz, Sz, y ) {
    i <- 1 : ( length(mz) - length(y) )
    m  <- mz[i]   + Sz[i,-i] %*% solve(Sz[-i,-i]) %*% (y-mz[-i])
    S  <- Sz[i,i] - Sz[i,-i] %*% solve(Sz[-i,-i]) %*% Sz[-i,i]
    return( list( m=m, S=S ) )
  }

  x1 <- c(0,4,4)    ; x2  <- c(3,0,3) ; y <- c(0.7,3.1,2.2) ; ny <- length(y)
  Vy <- 1.47        ; phi <- 1.9115
  x <- cbind(x1,x2) ; dx <- as.matrix( dist(x) ) ; Sy <- exp( -dx/phi ) * Vy

  mb   <- c(0,0,0) ; Sb <- diag(10 * 1.47,3) ; X <- cbind(1,x)
  b_y  <- GP.est( x, y, Sy, mb=mb, Sb=Sb, X=X )
  mb_y <- b_y$mb_y; Sb_y <- b_y$Sb_y
  # print( b_y )
  #  Mean: 1.2859283, 0.4238174, -0.2327936
  #  Covariance matrix 2.8847582 -0.46920017 -0.59063428
  #                   -0.4692002  0.14016302  0.06485419
  #                   -0.5906343  0.06485419  0.22413317

  x0   <- c(0,0) ; X0 <- c(1,x0)
  y0_y <- GP.pred( x0, x, y, Sy, phi, mb_y, Sb_y, X0=X0, X=X )
  m0_y <- y0_y$m0_y ; S0_y <- y0_y$S0_y     # 1.318117, 3.7413

# GP.pred( x[1,], x, y, Sy, phi, mb_y, Sb_y, X0=c(1,x[1,]), X=X )
# GP.pred( x[2,], x, y, Sy, phi, mb_y, Sb_y, X0=c(1,x[2,]), X=X )
# GP.pred( x[3,], x, y, Sy, phi, mb_y, Sb_y, X0=c(1,x[3,]), X=X )
# GP.pred( x0   , x, y, Sy, phi, mb_y, Sb_y, X0=X0        , X=X )
# GP.pred( x0   , x, y, Sy, phi, c(2,0,0), Sb*0, X0=X0    , X=X ) # Kriging (see GM)
# GP.pred( x0   , x, y, Sy, phi, c(0,0,0), Sb*0, X0=X0    , X=X ) # Kriging with 0-mean

  xy      <- as.geodata( cbind( x, y ) )
  xpred1  <- x0[1] ; xpred2 <- x0[2] ; xpred.1 <- expand.grid( xpred1, xpred2 )
  model.1 <- model.control( cov.m="exponential",
                            trend.d=~coords, trend.l=~coords )
  prior.1 <- prior.control( beta.prior="normal", beta=mb, beta.var.std=Sb/Vy,
                            sigmasq.prior="fixed", sigmasq=Vy,
                            phi.prior="fixed", phi=phi,
                            tausq.rel.prior="fixed",tausq.rel=0 )
  out.1   <- output.control(messages=FALSE)
  geoR.1  <- krige.bayes( xy, loc=xpred.1, mod=model.1, pr=prior.1, out=out.1 )

# with(geoR.1$posterior$beta$pars,{cat("Mean:",mean,"\nVariance:\n");print(var)})
# with(geoR.1$predictive,{cat("Mean:",mean,"\nVariance:\n");print(variance)})

# Geostatistical prediction after calibration on three data points marked by
# crosses. GP without nugget.
# Left: predictive means; right: predictive standard deviations.
  xpred.1 <- expand.grid(seq(-1,5,l=61), seq(-1,5,l=61))
  model.1 <- model.control( cov.m="exponential", trend.d=~coords, trend.l=~coords )
  prior.1 <- prior.control( beta.prior="normal", beta=mb, beta.var.std=Sb/Vy,
                             sigmasq.prior="fixed", sigmasq=Vy,
                             phi.prior="fixed", phi=phi,
                             tausq.rel.prior="fixed" , tausq.rel=0 )
  geoR.1  <- krige.bayes( xy, loc=xpred.1, mod=model.1, pr=prior.1, out=out.1 )
# We verify that the grid has the correct values (mean & var) at location (0,0):
  i_pred <- 621 # The index of pred.grid for which pred.grid[i_pred,] = c(0,0)
  # print( c(geoR.1$predictive$mean[i_pred],geoR.1$predictive$variance[i_pred])) # 1.318117,3.7413
# Then we create the regional maps
  par( mfrow=c(1,2) )
  image( geoR.1, val=geoR.1$predictive$mean, xlim=c(-2,6), ylim=c(-2,6),
        col=terrain.colors(21), main="\n\nMean", xlab="lon", ylab="lat" )
  points( xy, cex.max = 1, col = "black", add=T, pt.divide="equal", pch="x" )
  # contour( geoR.1, add=T, nlev=11 )
  legend.krige( val=geoR.1$predictive$mean,
                col=terrain.colors(21), x.leg=c(-1,5), y.leg=c(-2,-1.5) )
  image( geoR.1, val=sqrt(geoR.1$predictive$variance), xlim=c(-2,6), ylim=c(-2,6),
          col=terrain.colors(21), main="\n\nStandard Deviation", xlab="lon", ylab="lat" )
  points( xy, cex.max = 1, col = "black", add=T, pt.divide="equal", pch="x" )
  legend.krige( val=sqrt(geoR.1$predictive$variance),
                col=terrain.colors(21), x.leg=c(-1,5), y.leg=c(-2,-1.5) )

# iy0 <- which( xpred.1[,1] == x0[1] & xpred.1[,2]==x0[2] )
# iy1 <- which( xpred.1[,1] == x1[1] & xpred.1[,2]==x2[1] )
# iy2 <- which( xpred.1[,1] == x1[2] & xpred.1[,2]==x2[2] )
# iy3 <- which( xpred.1[,1] == x1[3] & xpred.1[,2]==x2[3] )
# print( c(geoR.1$predictive$mean[iy0], geoR.1$predictive$variance[iy0]) ) # 1.318117, 3.7413
# print( c(geoR.1$predictive$mean[iy1], geoR.1$predictive$variance[iy1]) ) # 0.7, 0
# print( c(geoR.1$predictive$mean[iy2], geoR.1$predictive$variance[iy2]) ) # 3.1, 0
# print( c(geoR.1$predictive$mean[iy3], geoR.1$predictive$variance[iy3]) ) # 2.2, 0

# Synug <- Sy + diag(diag(Sy)) * 0.5
# b_y   <- GP.est( x, y, Synug, mb=mb, Sb=Sb, X=X )
# mb_y  <- b_y$mb_y; Sb_y <- b_y$Sb_y
# GP.pred( x0   , x, y, Synug, phi, mb_y, Sb_y, X0=X0        , X=X )
# GP.pred( x[1,], x, y, Synug, phi, mb_y, Sb_y, X0=c(1,x[1,]), X=X )
# GP.pred( x[2,], x, y, Synug, phi, mb_y, Sb_y, X0=c(1,x[2,]), X=X )
# GP.pred( x[3,], x, y, Synug, phi, mb_y, Sb_y, X0=c(1,x[3,]), X=X )

  xpred.2 <- expand.grid(seq(-1,5,l=61), seq(-1,5,l=61))
  model.2 <- model.control( cov.m="exponential", trend.d=~coords, trend.l=~coords )
  prior.2 <- prior.control( beta.prior="normal", beta=mb, beta.var.std=Sb/Vy,
                            sigmasq.prior="fixed", sigmasq=Vy,
                            phi.prior="fixed", phi=phi,
                            tausq.rel.prior="fixed", tausq.rel=0.5 )

# Geostatistical prediction after calibration on three data points marked by
# crosses. GP with nugget.
# Left: predictive means; right: predictive standard deviations.
  geoR.2 <- krige.bayes( xy, loc=xpred.2, mod=model.2, pr=prior.2, out=out.1 )
  par( mfrow=c(1,2) )
  image( geoR.2, val=geoR.2$predictive$mean, xlim=c(-2,6), ylim=c(-2,6),
        col=terrain.colors(21), main="\n\nMean", xlab="lon", ylab="lat" )
  points( xy, cex.max = 1, col = "black", add=T, pt.divide="equal", pch="x" )
  # contour( geoR.2, add=T, nlev=11 )
  legend.krige( val=geoR.2$predictive$mean,
                col=terrain.colors(21), x.leg=c(-1,5), y.leg=c(-2,-1.5) )
  image( geoR.2, val=sqrt(geoR.2$predictive$variance), xlim=c(-2,6), ylim=c(-2,6),
          col=terrain.colors(21), main="\n\nStandard Deviation", xlab="lon", ylab="lat" )
  points( xy, cex.max = 1, col = "black", add=T, pt.divide="equal", pch="x" )
  legend.krige( val=sqrt(geoR.2$predictive$variance), scale.vals=c(0,1,2,3),
                col=terrain.colors(21), x.leg=c(-1,5), y.leg=c(-2,-1.5) )
  # print( c(geoR.2$predictive$mean[iy0], geoR.2$predictive$variance[iy0]) ) # 1.181519 , 5.02644
  # print( c(geoR.2$predictive$mean[iy1], geoR.2$predictive$variance[iy1]) ) # 0.6469166, 0
  # print( c(geoR.2$predictive$mean[iy2], geoR.2$predictive$variance[iy2]) ) # 3.038998 , 0
  # print( c(geoR.2$predictive$mean[iy3], geoR.2$predictive$variance[iy3]) ) # 2.255441 , 0

  # domain_phi <- phi * seq(0.1,3,0.1)
  # prior.3 <- prior.control( beta.prior="normal", beta=mb, beta.var.std=Sb/Vy,
  #   sigmasq.prior="sc.inv.chisq", sigmasq=Vy, df.sigmasq=2,
  #   phi.prior="uniform", phi.discrete=domain_phi,
  #   tausq.rel.prior="fixed", tausq.rel=0 )
  # geoR.3  <- krige.bayes( xy, loc=x0, mod=model.2, pr=prior.3, out=out.1 )
  # print(geoR.3$posterior$beta)    # Posterior(b)
  # print(geoR.3$posterior$sigmasq) # Posterior(sigmasq)
  # print(geoR.3$posterior$phi)     # Posterior(phi)
  # print(geoR.3$predictive)        # Prediction
  # print( sum( geoR.3$posterior$phi$phi.marginal$phi *
  #             geoR.3$posterior$phi$phi.marginal$expected ) ) # E[phi]
  # plot( domain_phi, geoR.3$posterior$phi$phi.marginal$expected )

# Upscaling error for three functions f(x). The mean of f(x) in the interval
# [-2,2] is calculated incorrectly as f applied to the mean of x (solid
# horizontal red line) and correctly as the average of f(x) (dashed horizontal
# blue line). Upscaling error f(E[x]) - E[f(x)] is negative for convex
# functions, zero for linear functions, positive for concave functions.
  x <- sort( runif(1e5,-2,2) )
  fConvex  <- function(x) { x^2 }
  fLinear  <- function(x) { x + 2 }
  fConcave <- function(x) { 4 - x^2 }
  yConvex <- fConvex(x) ; yLinear <- fLinear(x) ; yConcave <- fConcave(x) ;
  # cat( "m(x)=",mean(x),"; sd(x)=",sd(x),"; m(y)=",mean(y))
  par(mfrow=c(1,3), mar=c(2,3,2,1))
  # hist(x,main="")
  plot(x, yConvex , type="l", ylab="f(x)", lwd=1.5, bty="l")
  abline(h=fConvex(mean(x)),col="red") ; abline(h=mean(yConvex),lty=2,col="blue")
  legend( "top", title="f: CONVEX",
           legend=c("f(E[x])","E[f(x)]"),
           col=c("red","blue"), lty=c(1,2), lwd=c(1,1), cex=0.8 )
  plot(x, yLinear , type="l", ylab="", yaxt="n", lwd=1.5, bty="l")
  abline(h=fLinear(mean(x)),col="red") ; abline(h=mean(yLinear),lty=2,col="blue")
  legend( "top", title="f: LINEAR",
           legend=c("f(E[x])","E[f(x)]"),
           col=c("red","blue"), lty=c(1,2), lwd=c(1,1), cex=0.8 )
  plot(x, yConcave, type="l", ylab="", yaxt="n", lwd=1.5, bty="l")
  abline(h=fConcave(mean(x)),col="red") ; abline(h=mean(yConcave),lty=2,col="blue")
  legend( "bottom", title="f: CONCAVE",
           legend=c("f(E[x])","E[f(x)]"),
           col=c("red","blue"), lty=c(1,2), lwd=c(1,1), cex=0.8 )

# Exercise 1. Scaling error.
  x <- rnorm(1e5, 5, sqrt(10)) ; z <- x^3 ; meanz <- mean(z)