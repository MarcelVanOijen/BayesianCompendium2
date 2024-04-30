## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 22. Time Series and Data Assimilation

  GP.pred.0 <- function(x0,x,y,Sy,phi) {
    C0   <- Sy[1] * exp( -abs(x-x0)/phi )
    m0_y <- t(C0) %*% solve(Sy) %*% y
    S0_y <- Sy[1] - t(C0) %*% solve(Sy) %*% C0
    return( list( "m0_y"=m0_y, "S0_y"=S0_y ) ) }

  set.seed(13)
  Vy  <- 3 ; phi <- 1e1 ; nseries <- 1e3 ; n <- 30
  Y       <- NULL
  for(i in 1:nseries) {
    x   <- 0 ; y  <- rnorm(1,0,sqrt(Vy))
    for(xi in 1:n) {
      dx <- as.matrix( dist(x) ) ; Sy <- exp( -dx/phi ) * Vy
      py <- GP.pred.0(xi,x,y,Sy,phi)
      yi <- rnorm(1, mean=py$m0_y, sd=sqrt(py$S0_y))
      x  <- c(x,xi) ; y <- c(y,yi)
    } ; Y <- cbind( Y, y )
  } ; mY <- rowMeans(Y) ; sY <- apply(Y, 1, sd)

  set.seed(13)
  alpha <- exp(-1/phi) ; Ve <- Vy * (1-alpha^2)
  Z <- NULL
  for(i in 1:nseries) {
    z <- rnorm( 1, 0, sqrt(Ve/(1-alpha^2)) )
    for(t in 1:n){
      z1 <- rnorm(1, mean=alpha*tail(z,1), sd=sqrt(Ve))
      z  <- c(z,z1)
    } ; Z <- cbind( Z, z )
  } ; mZ <- rowMeans(Z) ; sZ <- apply(Z, 1, sd)

# Left: sampling from a Gaussian Process. Black lines are mean plus-minus one
# standard deviation of 1000 sampled functions for each value of x. Coloured
# lines show three sampled functions. Right: as on the left but from an AR(1)
# process rather than a GP. The same random seed was used for both panels.
  par( mfrow=c(1,2), mar=c(5,4,4,2) )
  x <- 0:n ; yrange <- 2*range(mY-sY,mY+sY)
  plot  (x, mY           , type="l", lwd=3, ylim=yrange,
         ylab="y", main="GP")
  points(x, mY-sY        , type="l", lwd=3, lty=2 )
  points(x, mY+sY        , type="l", lwd=3, lty=2 )
  points(x, Y[,1        ], type="l", col="red"   )
  points(x, Y[,nseries/2], type="l", col="blue"  )
  points(x, Y[,nseries  ], type="l", col="green" )
  plot  (x, mZ           , type="l", lwd=3, ylim=yrange,
         xlab="t", ylab="z", main="AR(1)")
  points(x, mZ-sZ        , type="l", lwd=3, lty=2 )
  points(x, mZ+sZ        , type="l", lwd=3, lty=2 )
  points(x, Z[,1        ], type="l", col="red"   )
  points(x, Z[,nseries/2], type="l", col="blue"  )
  points(x, Z[,nseries  ], type="l", col="green" )

  set.seed(13)
  n     <- 30 ; y <- runif(n,0.8,1) ; Vy <- 3
  zp    <- Vp <- za <- Va <- rep(NA,n)
  zp[1] <- 0 ; Vp[1] <- Vf <- 1 ; alpha <- 0.9
  for(i in 1:n) {
    Ki      <- Vp[i] / ( Vp[i] + Vy )
    za[i]   <- zp[i] + Ki * ( y[i] - zp[i] )
    Va[i]   <- ( 1 - Ki ) * Vp[i] ; if(i==n) break
    zp[i+1] <- alpha   * za[i]
    Vp[i+1] <- alpha^2 * Va[i] + Vf }

# Data assimilation using the KF. The top panel shows the observations and the
# alternating state estimates from model prediction (mp) and data assimilation
# (ma). The bottom panel shows the corresponding uncertainties as variances.
# Observational uncertainty is not shown but has a constant variance equal to 3.
  tp  <- (1:n) ; ta <- tp + 1e-3
  tpa <- sort( c(tp,ta) ) ; zpa <- c( rbind(zp,za) ) ; Vpa <- c( rbind(Vp,Va) )
  par( mfrow=c(2,1), mar=c(1,2,1,0) )
  plot  ( 1:n, zp, ylim=c(0,1), xlab="time",ylab="", col="red", xaxt="n" )
  points( 1:n, za, col="black", pch=16)
  points( 1:n,  y, col="blue", pch=4)
  points( tpa, zpa, type="l" )
  legend( "bottomright", title=paste0("Mean"),
          legend=c("data","ma","mp"),
          col=c("blue","black","red"), pch=c(4,16,1) )
  par( mar=c(2,2,0,0) )
  plot  ( 1:n, Vp, ylim=c(0,2), xlab="time",ylab="", col="red" )
  points( 1:n, Va, col="black", pch=16)
  points( tpa, Vpa, type="l" )
  legend( "bottomright", title=paste0("Variance"),
          legend=c("Vp","Va"),
          col=c("red","black"), pch=c(1,16) )

  # Vy <- rep(1e8,n)
  # alpha <- exp(-1/10) ; Vf <- 3 * (1-alpha^2) ; zp[1] <- 0 ; Vp[1] <- 3
  # for(i in 1:n) {
  #   Ki      <- Vp[i] / ( Vp[i] + Vy[i] )
  #   za[i]   <- zp[i] + Ki * ( y[i] - zp[i] )
  #   Va[i]   <- ( 1 - Ki ) * Vp[i] ; if(i==n) break
  #   zp[i+1] <- alpha   * za[i]
  #   Vp[i+1] <- alpha^2 * Va[i] + Vf }
  # print(zp) ; print(Vp) ; print(za) ; print(Va)

  # set.seed(13)
  # n <- 30 ; y <- runif(n,0.8,1) ; Vy <- rep(3,n)
  # zp <- Sp <- za <- Sa <- vector("list",n)
  # zp[[1]]  <- matrix(0) ; Sp[[1]] <- Vf <- matrix(1) ; F <- matrix(0.9)
  # H        <- matrix(1)
  # for(i in 1:n) {
  #   Ki        <- Sp[[i]] %*% t(H) %*% solve( H %*% Sp[[i]] %*% t(H) + Vy[i] )
  #   za[[i]]   <- zp[[i]] + Ki %*% ( y[i] - H %*% zp[[i]] )
  #   Sa[[i]]   <- ( diag(1) - Ki %*% H ) %*% Sp[[i]] ; if(i==n) break
  #   zp[[i+1]] <- F %*% za[[i]]
  #   Sp[[i+1]] <- F %*% Sa[[i]] %*% t(F) + Vf }
  # print(unlist( zp )) ; print(unlist( Sp ))
  # print(unlist( za )) ; print(unlist( Sa ))

  set.seed(13)
  n  <- 30 ; y <- runif(n,0.8,1) ; Vy <- rep(3,n)

  zp      <- Sp <- za <- Sa <- vector("list",n)
  zp[[1]] <- c(0,0) ; Sp[[1]] <- Sf <- diag(1,2)
  F       <- diag(0.9,2)
  H       <- t( c(0.5,0.5) )
  for(i in 1:n) {
    Ki        <- Sp[[i]] %*% t(H) %*% solve( H %*% Sp[[i]] %*% t(H) + Vy[i] )
    za[[i]]   <- zp[[i]] + Ki %*% ( y[i] - H %*% zp[[i]] )
    Sa[[i]]   <- ( diag(2) - Ki %*% H ) %*% Sp[[i]] ; if(i==n) break
    zp[[i+1]] <- F %*% za[[i]]
    Sp[[i+1]] <- F %*% Sa[[i]] %*% t(F) + Sf }

  # print(unlist( zp )) ; print(unlist( Sp ))
  # print(unlist( za )) ; print(unlist( Sa ))
  # print( sum(Sp[[30]])/4 )
  # print( sum(Sa[[30]])/4 )

# Data assimilation using a bivariate KF. The top panel shows the observations
# and the alternating mean state estimates from model prediction (mp) and data
# assimilation (ma). The bottom panel shows the corresponding uncertainties.
# Observational uncertainty is not shown but has a constant variance equal to 3.
  tp   <- (1:n) ; ta <- tp + 1e-3 ; tpa <- sort( c(tp,ta) )
  zp.u <- unlist(zp) ; Sp.u <- unlist(Sp) ; za.u <- unlist(za) ; Sa.u <- unlist(Sa)
  iz1  <- seq(1, 59,by=2) ; iz2  <- iz1+1
  iV1  <- seq(1,117,by=4) ; iC12 <- iV1+1 ; iV2 <- iV1+3
  mp   <- 0.5 * (zp.u[iz1] + zp.u[iz2])
  ma   <- 0.5 * (za.u[iz1] + za.u[iz2])
  Vp   <- (Sp.u[iV1] + Sp.u[iV2] + 2 * Sp.u[iC12]) * 0.5^2
  Va   <- (Sa.u[iV1] + Sa.u[iV2] + 2 * Sa.u[iC12]) * 0.5^2
  mpa  <- c( rbind(mp,ma) ) ; Vpa <- c( rbind(Vp,Va) )
  par( mfrow=c(2,1), mar=c(1,2,1,0) )
  plot  ( tp ,  mp, ylim=c(0,1), xlab="time",ylab="", col="red", xaxt="n" )
  points( ta ,  ma, col="black", pch=16)
  points( tp ,   y, col="blue", pch=4)
  points( tpa, mpa, type="l" )
  legend( "bottomright", title=paste0("Mean"),
          legend=c("data","ma","mp"),
          col=c("blue","black","red"), pch=c(4,16,1) )
  par( mar=c(2,2,0,0) )
  plot  ( tp , Vp , ylim=c(0,2), xlab="time",ylab="", col="red" )
  points( ta , Va , col="black", pch=16)
  points( tpa, Vpa, type="l" )
  legend( "bottomright", title=paste0("Variance"),
          legend=c("Vp","Va"),
          col=c("red","black"), pch=c(1,16) )

  KalmanGain_SVD <- function(Sp,H,Vy) {
    Z     <- H %*% Sp %*% t(H) + Vy ; svd_Z <- svd(Z)
    inv_Z <- solve( t(svd_Z$v) ) %*% solve( diag(svd_Z$d) ) %*% solve( svd_Z$u )
    return( Sp %*% t(H) %*% inv_Z ) }
  