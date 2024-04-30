## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 7. Sampling from the Posterior Distribution by MCMC

  library(mvtnorm)

  dnorm( 38:40, m=0, sd=1 ) ; dnorm( 38:40, m=0, sd=1, log=T)

  z <- rep( 3, 200 ) ; prod( dnorm(z) ) ; sum ( dnorm(z,log=T) )

  x  <- c(10,20,30) ; y  <- c(6.09,8.81,10.66)
  ny <- length(y)   ; Sy <- diag(3,ny)

  f  <- function( x, b ) { cbind(1,x)%*%b }
  mb <- c(0,0) ; nb <- length(mb) ; Sb <- diag(1e4,nb)
  logPrior <- function( b, mb.=mb, Sb.=Sb )
    { dmvnorm( b, mb., Sb., log=T ) }
  logLik   <- function( b, y.=y, Sy.=Sy, x.=x, f.=f ) {
    dmvnorm( y., mean=f.(x, b), sigma=Sy., log=T ) }
  logPost <- function( b, mb.=mb, Sb.=Sb, y.=y, Sy.=Sy, x.=x, f.=f ) {
    logPrior(b, mb., Sb.) + logLik(b, y., Sy., x., f.) }

  MetropolisLogPost <- function( f.=f, logp=logPost, b0=mb,
                                 SProp=Sb/100, ni=1e4, ... ) {
    bChain <- matrix( NA, nrow=ni, ncol=length(b0) ) ; bChain[1,] <- b0
    for( i in 2:ni ) {
       b1  <- bChain[i-1,] + as.numeric( rmvnorm(1,sigma=SProp) )
       if( log(runif(1)) < logp(b1,f.=f.,...) - logp(b0,f.=f.,...) )
            bChain[i,] <- b0 <- b1
       else bChain[i,] <- b0 }
    return(bChain) }

  bChain1 <- MetropolisLogPost( f )

# Failed Metropolis sampling of the straight-line model using the default settings for the proposal distribution.
  par( mfrow=c(2,nb), mar=c(2,2,3,2) )
  pNames <- c( "Intercept", "Slope" ) ; pMeans <- colMeans(bChain1)
  for(i in 1:nb){ plot( bChain1[,i], xlab="", ylab="", pch=".",
                        main=paste(pNames[i],"\nmean=",signif(pMeans[i],3) ) ) }
  for(i in 1:nb){ hist( bChain1[,i], xlab="", ylab="", main=pNames[i] ) }

  bChain2  <- MetropolisLogPost( f, SProp=diag(c(0.1,0.01)) )

# Metropolis sampling of the straight-line model with small proposal variances.
  par( mfrow=c(2,nb), mar=c(2,2,3,2) )
  pNames <- c( "Intercept", "Slope" ) ; pMeans <- colMeans(bChain2)
  for(i in 1:nb){ plot( bChain2[,i], xlab="", ylab="", pch=".",
                        main=paste(pNames[i],"\nmean=",signif(pMeans[i],3) ) ) }
  for(i in 1:nb){ hist( bChain2[,i], xlab="", ylab="", main=pNames[i] ) }

  set.seed(13)

  x       <- c(10,20,30) ; yq <- y + 0.1*x^2
  nyq     <- length(yq)  ; Syq <- diag(30,nyq)
  fq      <- function( x, b ) { cbind(1,x,x^2)%*%b }
  nbq     <- 3 ; mbq <- rep(0,nbq) ; Sbq <- diag(nbq) * 1e4
  bqChain <- MetropolisLogPost( fq, b0=c(0,0,0),
    SProp=diag(c(1,0.1,0.01)), mb.=mbq, Sb.=Sbq, y.=yq )

# Metropolis sampling of the 3-parameter quadratic model.
  par( mfrow=c(2,nbq), mar=c(2,2,3,2) )
  pNames <- c( "Intercept", "Slope", "Slope2" ) ; pMeans <- colMeans(bqChain)
  for(i in 1:nbq){ plot( bqChain[,i], xlab="", ylab="", pch=".",
                        main=paste(pNames[i],"\nmean=",signif(pMeans[i],3) ) ) }
  for(i in 1:nbq){ hist( bqChain[,i], xlab="", ylab="", main=pNames[i] ) }

# Scatter plots for the two MCMCs of the linear model."}
  par(mfrow=c(1,2))
  plot(bChain1, type="b", xlab="b1", ylab="b2", pch=".", main="Default")
  plot(bChain2, type="b", xlab="b1", ylab="b2", pch=".", main="Small prop. var.")

# Twenty predictions from the posterior distribution for the linear model."}
  par( mfrow=c(1,1) )
  plot( x, f(x,bChain2[100,]), type="l", ylim=c(0,15), ylab="y and f(x)" )
  points( cbind(x,y), col='blue', cex=2, lwd=3 )
  isubsample <- round( seq(1e4/2,1e4,length.out=20) )
  for(i in isubsample) { lines( x, f(x,bChain2[i,]), lty="dashed" ) }
  
# EXERCISE 1.
  f2      <- function( x, b ) { cbind(1,x)%*%b }
  mb2     <- c(4,0.2) ; Sb2 <- diag(1,2)
  SProp2  <- diag(0.1,2)
  bChain2 <- MetropolisLogPost( f2, b0=mb2, SProp=SProp2, mb=mb2, Sb=Sb2 )
  
  f3      <- function( x, b ) { cbind(1,x,x)%*%b }
  mb3     <- c(4,0.1,0.1) ; Sb3 <- diag(c(1,1/2,1/2))
  SProp3  <- diag(c(0.1,0.1/2,0.1/2))
  bChain3 <- MetropolisLogPost( f3, b0=mb3, SProp=SProp3, mb=mb3, Sb=Sb3 )

# EXERCISE 3: Try matching panels A-F to MCMCs 1-6 described in the text."}
  mb <- c(4,0.2) ; Sb <- diag(1,2) ; Sy <- diag(3,3) ; SProp <- diag(0.1,2) ; ni <- 5e3
  bChains      <- vector( "list", 6 )
  bChains[[1]] <- MetropolisLogPost( f, Sb=Sb, Sy=Sy, SProp=SProp, ni=ni )
  bChains[[3]] <- MetropolisLogPost( f, Sb=Sb, Sy=Sy*10, SProp=SProp, ni=ni )
  bChains[[5]] <- MetropolisLogPost( f, Sb=Sb/c(100,1), Sy=Sy, SProp=SProp, ni=ni )
  bChains[[6]] <- MetropolisLogPost( f, Sb=Sb/5, Sy=Sy, SProp=SProp, ni=ni )
  bChains[[4]] <- MetropolisLogPost( f, Sb=Sb, Sy=Sy, SProp=SProp*50, ni=ni )
  bChains[[2]] <- MetropolisLogPost( f, Sb=Sb, Sy=Sy, SProp=SProp, ni=ni/100 )
  par( mfrow=c(2,3), mar=c(2,2,3,2) )
  for( im in 1:6 ) { plot( bChains[[im]], type="b", xlab="b1", ylab="b2", pch=".",
    main=paste("Chain",LETTERS[im]), xlim=c(1,7), ylim=c(-0.1,0.6) ) }
