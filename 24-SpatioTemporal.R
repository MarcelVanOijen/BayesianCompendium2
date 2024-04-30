## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 24. Spatio-Temporal Modelling and Adaptive Sampling

  library(raster)

  coords   <- as.matrix( expand.grid(1:7,1:7) ) ; nz  <- dim(coords)[1]
  tausq    <- 10 ; Vy <- diag( tausq, 1 ) ; ny <- 8 ; y <- 0
  mb       <-  0 ; Va <- 80 ; rho <- 0.75
  Rspatial <- rho^as.matrix( dist(coords) ) 
  zp       <- vector("list",ny)   ; Sp      <- vector("list",ny)
  za       <- vector("list",ny+1) ; Sa      <- vector("list",ny+1)
  za[[1]]  <- rep( mb, nz )       ; Sa[[1]] <- Va * Rspatial
  APV      <- rep(NA,ny+1)        ; APV[1]  <- mean( diag( Sa[[1]] ) )
  alpha    <- 0.9                 ; M <- diag(alpha,nz)
  Vf       <- 20                  ; Q <- Vf * Rspatial
  iSensor  <- rep(NA,ny)

  KalmanGain <- function(Sp,H,Vy){ Sp %*% t(H) %*% solve(H%*%Sp%*%t(H) + Vy) }
  for (i in 1:ny) {
    zp[[i]] <- M %*% za[[i]] ; Sp[[i]] <- M %*% Sa[[i]] %*% t(M) + Q
    r_APV <- rep(NA,nz)
    for (r in 1:nz) {
      r_Hi     <- matrix(0,nrow=1,ncol=nz) ; r_Hi[r] <- 1
      r_Ki     <- KalmanGain( Sp[[i]], r_Hi, Vy )
      r_Sa     <- ( diag(nz) - r_Ki %*% r_Hi ) %*% Sp[[i]]
      r_APV[r] <- mean( diag(r_Sa) ) }
    iSensor[i] <- match( min(r_APV), r_APV )
    Hi         <- matrix(0,nrow=1,ncol=nz) ; Hi[iSensor[i]] <- 1
    Ki         <- KalmanGain( Sp[[i]], Hi, Vy )
    za[[i+1]]  <- zp[[i]] + Ki %*% ( y - Hi %*% zp[[i]] )
    Sa[[i+1]]  <- ( diag(nz) - Ki %*% Hi ) %*% Sp[[i]]
    APV[i+1]   <- mean( diag( Sa[[i+1]] ) ) }

  par( mfrow=c(2,5), mar=c(0,0,1,0) )
  for(i in 1:8) {
    plot( rasterFromXYZ(cbind(coords, diag(Sa[[i]]) )),
          main=paste0("t = ",i-1), breaks=seq(20,90,length.out=8),
          col=gray.colors(8,rev=T), axes=FALSE, box=FALSE, legend=FALSE )
    if(i>1) { points(coords[iSensor[i-1],1],coords[iSensor[i-1],2],pch="x")}
  }
  plot( rasterFromXYZ(cbind(coords, diag(Sa[[9]]) )),
        main=paste0("t = ",8), breaks=seq(20,90,length.out=8),
        col=gray.colors(8,rev=T), axes=FALSE, box=FALSE )
  points(coords[iSensor[9-1],1],coords[iSensor[9-1],2],pch="x")
