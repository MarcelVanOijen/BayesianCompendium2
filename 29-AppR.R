## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Appendix D: R

  x    <- 2
  y    <- c(1,2)
  yrow <- matrix( c(1,2), nrow=1) ; ycol <- matrix( c(1,2), ncol=1)
  Z    <- matrix( 0:3, nrow=2)

  print( Z ) ; print( length(Z) ) ; print( dim(Z) )

  xZ <- x * Z
  Zy <- Z %*% y ; yZ <- y %*% Z ; Zyc <- Z %*% ycol ; yrZ <- yrow %*% Z

  I5a <- diag( 5 ) ; I5b <- diag( rep(1,5) ) ; X <- diag( 1:3 )

  invZ <- solve(Z) ; trZ <- t(Z)

  mylist <- list( distr="Bivariate Gaussian", m=c(0,0), S=diag(2) )
  print( mylist$distr ) ; print( mylist[ 2:3 ] )

  myf <- function(x){ x^2 }

  curve( myf )
  plot( 1:10, myf(1:10), main="My function", xlab="x", ylab="y" )

  x  <- rnorm( 1e4, mean=1, sd=0.3 )
  dx <- dnorm( x  , mean=1, sd=0.3 )

  hist ( x )
  curve( dnorm,-2,2 )
  curve( dbeta(x,4,2) )

  barplot( dbinom( 0: 1, size= 1, p=0.7 ), names.arg=0: 1, main="Br[x]" )
  barplot( dbinom( 0:10, size=10, p=0.7 ), names.arg=0:10, main="Bi[x]" )
  barplot( dpois ( 0:10, lambda=1 )      , names.arg=0:10, main="Pn[x]" )

  sample <- rmvnorm( 100, sigma=diag(1,2) ) # Gives error message
  # install.packages( 'mvtnorm' )           # Do this only once
  library( mvtnorm )                        # Do this each new R-session
  sample <- rmvnorm( 100, sigma=diag(1,2) ) # Now it works
  plot( sample )
  