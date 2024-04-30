## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 4. Assigning a Likelihood Function

  library(mvtnorm)

  x <- 1:3 ; y <- c(10,1,5) ; ny <- length(y)
  X    <- cbind( 1, x ) # Design matrix

  Vy1 <- rep(1,ny) ; Sy1 <- diag(Vy1) # Var constant
  Vy2 <- y/max(y)  ; Sy2 <- diag(Vy2) # Var prop. to y
  Vy3 <- (y/10)^2  ; Sy3 <- diag(Vy3) # Var prop. to y^2 (= CV constant)

  Sb_y1 <- solve(    t(X) %*% solve(Sy1) %*% X ) # GLS variance
  mb_y1 <- Sb_y1 %*% t(X) %*% solve(Sy1) %*% y   # GLS estimator
  
  Sb_y2 <- solve(    t(X) %*% solve(Sy2) %*% X ) # GLS variance
  mb_y2 <- Sb_y2 %*% t(X) %*% solve(Sy2) %*% y   # GLS estimator
  
  Sb_y3 <- solve(    t(X) %*% solve(Sy3) %*% X ) # GLS variance
  mb_y3 <- Sb_y3 %*% t(X) %*% solve(Sy3) %*% y   # GLS estimator

# Posterior mean straight lines derived by Bayesian calibration using three
# different likelihood functions.  Measurement variance is constant (Case 1),
# proportional to $y$ (Case 2), proportional to $y^2$ (Case 3).
  par( mfrow=c(1,1), mar=c(4,4,4,4) )
  plot(x,y,xlim=c(0.8,3.2),ylim=c(0,11))
  points( x, X %*% mb_y1, type="l", col="black")
  points( x, X %*% mb_y2, type="l", col="red")
  points( x, X %*% mb_y3, type="l", col="blue")
  legend( "topright", legend=c("Case 1","Case 2","Case 3"),
          col=c("black","red","blue"), lty=1 )
  