## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 6. Markov Chain Monte Carlo Sampling (MCMC)

  library(mvtnorm)

  set.seed(13)

  ni <- 10 ; b <- rep(NA,ni) ; b[1] <- 0 ; pjump <- 0.5
  for(i in 2:ni){ jump <- (runif(1)<pjump) ; b[i] <- (b[i-1]+jump)%%2 }

  set.seed(13)
  ni2 <- 1e3 ; b2 <- rep(NA,ni2) ; b2[1] <- 0 ; pjump <- 0.5
  for(i in 2:ni2){ jump <- (runif(1)<pjump) ; b2[i] <- (b2[i-1]+jump)%%2 }

# Sampling from a Bernoulli distribution using MCMC. Top panels: MCMC with 10 iterations; bottom panels: 1000 iterations. Left: trace plots; right: chain histograms.
  par(mfrow=c(2,2),mar=c(2,2,2,2))
  plot(b, xlab="Iteration (i)",ylab="b[i]", main="Trace plot")
  hist(b,10, xaxt='n') ; axis(side=1, at=c(0.05,0.95), labels=0:1)
  plot(b2, xlab="Iteration (i)",ylab="", main="")
  hist(b2,10, xaxt='n', main="") ; axis(side=1, at=c(0.05,0.95), labels=0:1)

# Sampling from a Bernoulli distribution using a poorly mixing (jump probability
# 0.01) MCMC with 1000 iterations. Left: trace plot; right: chain histogram.
  set.seed(1)
  par( mfrow=c(1,2), mar=c(2,2,2,2) )
  ni3 <- 1e3 ; b3 <- rep(NA,ni3) ; b3[1] <- 0 ; pjump <- 0.01
  for(i in 2:ni3){ jump <- (runif(1)<pjump) ; b3[i] <- (b3[i-1]+jump)%%2 }
  plot( b3, xlab="Iteration (i)",ylab="", main="")
  hist( b3, 10, xaxt='n', main="" )
  axis( side=1, at=c(0.05,0.95), labels=0:1 )

  Metropolis <- function( p, b0=mb, SProp=Sb/100, ni=1e4 ){
    bChain <- matrix( NA, nrow=ni, ncol=length(b0) ) ; bChain[1,] <- b0
    for( i in 2:ni ){
       b1 <- bChain[i-1,] + as.numeric( rmvnorm(1,sigma=SProp) )
       if( runif(1) < p(b1)/p(b0) ) bChain[i,] <- b0 <- b1
       else                         bChain[i,] <- b0 }
    return(bChain) }

  fA <- function(b){ exp(-b^2) }
  bA <- Metropolis( fA, b0=0, SProp=diag(1) )
  mA <- round( mean(bA), 2 ) ; VA <- round( var(bA), 2 )
# Sampling from an unknown distribution using the Metropolis algorithm, where
# all we know is that p[b] is proportional to $\\exp (-b^2)$.
  par( mfrow=c(1,2), mar=c(2,2,2,2) )
  plot( bA, xlab="Iteration (i)", ylab="", main="Trace plot" )
  hist( bA, main=paste0( "m=", mA, ", V=",VA) )

  fB <- function(b){ exp(-b[1]^2 - b[2]^2) }
  bB <- Metropolis( fB, b0=c(0,0), SProp=diag(2)  )
# Sampling from an unknown bivariate distribution using the Metropolis algorithm (see text).
  par( mfrow=c(1,3) )
  plot( bB, asp=1, type="b", main="2D trace plot" )
  hist( bB[,1], main="b1" ) ; hist( bB[,2], main="b2" )

  fC <- function(b){ b[1]^2 + b[2]^2 < 1 || (b[1]-2)^2 + (b[2]-0.5)^2 < 0.5 }
  bC <- Metropolis( fC, b0=c(0,0), SProp=diag(10,2) )
# Sampling from an unknown bimodal, bivariate distribution using the Metropolis
# algorithm (see text).
  par( mfrow=c(1,3) )
  plot( bC, asp=1, type="p", main="2D trace plot" )
  hist( bC[,1], main="b1") ; hist( bC[,2], main="b2" )

