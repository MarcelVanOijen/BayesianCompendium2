## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 21. Machine Learning

  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(mvtnorm)
  library(rsvg)

# Family tree of machine learning methods.
  DAGML <- grViz( "digraph{ graph[]
    node[ shape=box ]
      A [label='@@1'] ; B1[label='@@2'] ; B2[label='@@3']
      C1[label='@@4'] ; C2[label='@@5'] ; C3[label='@@6'] ; C4[label='@@7']
    edge[]
      A -> B1 ; A -> B2 ; B1 -> C1 ; B1 -> C2 ; B2 -> C3 ; B2 -> C4 }
    [1]: 'Machine learning'
    [2]: 'Supervised'
    [3]: 'Unsupervised'
    [4]: 'Regression'
    [5]: 'Classification'
    [6]: 'Clustering'
    [7]: paste0( 'Dimensionality', '\\nreduction' )
  ")
  export_svg(DAGML) %>% charToRaw() %>% rsvg() %>% png::writePNG("DAGML.png")

# A 2-3-2 neural network. Each neuron is connected to all inputs and to all
# outputs and there are no loops, making this a standard feedforward neural
# network.
  DAGNN <- grViz("digraph{ graph[rankdir=LR, ranksep=1, nodesep=0.5]
    subgraph cluster_0 {
      graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='OUTPUTS']
      node[fillcolor=White, margin=0.1, shape=circle, style=filled]
      y1[label = '@@6']
      y2[label = '@@7']
    }
    subgraph cluster_1 {
      graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='HIDDEN LAYER']
      node[fillcolor=White, margin=0.1, shape=circle, style=filled]
      h1[label = '@@4']
      h2[label = '@@5']
      h3[label = '@@3']
    }
    subgraph cluster_2 {
      graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='INPUTS']
      node[fillcolor=White, margin=0.1, shape=circle, style=filled]
      x1[label = '@@1']
      x2[label = '@@2']
    }
    edge[arrowhead=vee, arrowsize=1.25]
      x1 -> h1 ; x1 -> h2 ; x1 -> h3
      x2 -> h1 ; x2 -> h2 ; x2 -> h3
      h1 -> y1 ; h1 -> y2
      h2 -> y1 ; h2 -> y2
      h3 -> y1 ; h3 -> y2 }
    [1]: paste0('x&#x2081;')
    [2]: paste0('x&#x2082;')
    [3]: paste0('h&#x2081;')
    [4]: paste0('h&#x2082;')
    [5]: paste0('h&#x2083;')
    [6]: paste0('y&#x2081;')
    [7]: paste0('y&#x2082;')
  ")
  export_svg(DAGNN) %>% charToRaw() %>% rsvg() %>% png::writePNG("DAGNN.png")

  fNN.1n1 <- function( x, b1, b2, w1, w2, fh=tanh, fy=identity ) {
    X      <- cbind( 1, x )
    theta1 <- rbind( b1, w1 )    ; theta2 <- c( b2, w2 )
    h      <- fh( X %*% theta1 ) ; H      <- cbind( 1, h )
    y      <- fy( H %*% theta2 )
    return( y ) }

  fNN.1nn1 <- function( x, b1, b2, b3, w1, w2, w3, fh=tanh, fy=identity ) {
    n1     <- length( b1 ) ; n2 <- length( b2 )
    X      <- cbind( 1, x )
    theta1 <- rbind( b1, w1 )
    theta2 <- matrix( c( b2, w2 ), nrow=n1+1, byrow=T )
    theta3 <- c( b3, w3 )
    h1     <- fh( X  %*% theta1 ) ; H1 <- cbind( 1, h1 )
    h2     <- fh( H1 %*% theta2 ) ; H2 <- cbind( 1, h2 )
    y      <- fy( H2 %*% theta3 )
    return( y ) }

# Sampling from different neural network architectures, as specified in the
# headers, with the number of parameters between brackets. The number of hidden
# layers was one (top row) or two (bottom). Black lines: mean of 500 network
# functions sampled from the prior $\\pm$ standard deviation. Priors were
# $N[\\mu=0,\\sigma=3]$ for all parameters. Coloured lines: three of the sampled
# network functions.
  nx    <- 500 ; x   <- seq(-5,5,length.out=nx)
  nout  <- 1e3 ; out <- matrix(NA,nrow=nx,ncol=nout)
  sd_w1 <- sd_b1 <- sd_w2 <- sd_b2 <- sd_w3 <- sd_b3 <- 3
  set.seed(1)
  par(mfrow=c(2,3), mar=c(2,3,2,1))
  # 1-4-1
  n1 <- 4
  for( j in 1:nout ) {
    w1 <- rnorm( n1, 0, sd_w1 )
    b1 <- rnorm( n1, 0, sd_b1 )
    w2 <- rnorm( n1, 0, sd_w2 )
    b2 <- rnorm( 1 , 0, sd_b2 )
    out[,j] <- fNN.1n1( x, b1, b2, w1, w2 )
    # out[,j] <- fNN.1n1( x, b1, b2, w1, w2, fh=plogis )
  }
  m <- rowMeans(out) ; s <- apply(out, 1, sd)
  plot  (x, m           , type="l", lwd=3,
                          main=expression("1-4-1 (n"[p]*" = 13)"),
                          xaxt='n', xlab="", ylab="", ylim=c(-15,15) )
  points(x, m-s         , type="l", lty=2 )
  points(x, m+s         , type="l", lty=2 )
  points(x, out[,1     ], type="l", col="red"   )
  points(x, out[,nout/2], type="l", col="blue"  )
  points(x, out[,nout  ], type="l", col="green" )
  # 1-20-1
  n1 <- 20
  for( j in 1:nout ) {
    w1 <- rnorm( n1, 0, sd_w1 )
    b1 <- rnorm( n1, 0, sd_b1 )
    w2 <- rnorm( n1, 0, sd_w2 )
    b2 <- rnorm( 1 , 0, sd_b2 )
    out[,j] <- fNN.1n1( x, b1, b2, w1, w2 )
  }
  m <- rowMeans(out) ; s <- apply(out, 1, sd)
  plot  (x, m           , type="l", lwd=3,
                          main=expression("1-20-1 (n"[p]*" = 61)"),
                          xaxt='n', xlab="", ylab="", ylim=c(-30,30) )
  points(x, m-s         , type="l", lty=2 )
  points(x, m+s         , type="l", lty=2 )
  points(x, out[,1     ], type="l", col="red"   )
  points(x, out[,nout/2], type="l", col="blue"  )
  points(x, out[,nout  ], type="l", col="green" )
  plot.new()
  # 1-2-2-1
  n1 <- n2 <- 2
  for( j in 1:nout ) {
    w1 <- rnorm( n1   , 0, sd_w1 )
    b1 <- rnorm( n1   , 0, sd_b1 )
    w2 <- rnorm( n1*n2, 0, sd_w2 )
    b2 <- rnorm(    n2, 0, sd_b2 )
    w3 <- rnorm(    n2, 0, sd_w3 )
    b3 <- rnorm(     1, 0, sd_b3 )
    out[,j] <- fNN.1nn1( x, b1, b2, b3, w1, w2, w3 )
  }
  m <- rowMeans(out) ; s <- apply(out, 1, sd)
  plot  (x, m           , type="l", lwd=3,
                          main=expression("1-2-2-1 (n"[p]*" = 13)"),
                          xlab="", ylab="", ylim=c(-15,15) )
  points(x, m-s         , type="l", lty=2 )
  points(x, m+s         , type="l", lty=2 )
  points(x, out[,1     ], type="l", col="red"   )
  points(x, out[,nout/2], type="l", col="blue"  )
  points(x, out[,nout  ], type="l", col="green" )
  # 1-6-6-1
  n1 <- n2 <- 6
  for( j in 1:nout ) {
    w1 <- rnorm( n1   , 0, sd_w1 )
    b1 <- rnorm( n1   , 0, sd_b1 )
    w2 <- rnorm( n1*n2, 0, sd_w2 )
    b2 <- rnorm(    n2, 0, sd_b2 )
    w3 <- rnorm(    n2, 0, sd_w3 )
    b3 <- rnorm(     1, 0, sd_b3 )
    out[,j] <- fNN.1nn1( x, b1, b2, b3, w1, w2, w3 )
  }
  m <- rowMeans(out) ; s <- apply(out, 1, sd)
  plot  (x, m           , type="l", lwd=3,
                          main=expression("1-6-6-1 (n"[p]*" = 61)"),
                          xlab="", ylab="", ylim=c(-30,30) )
  points(x, m-s         , type="l", lty=2 )
  points(x, m+s         , type="l", lty=2 )
  points(x, out[,1     ], type="l", col="red"   )
  points(x, out[,nout/2], type="l", col="blue"  )
  points(x, out[,nout  ], type="l", col="green" )
  # 1-10-10-1
  n1 <- n2 <- 10
  for( j in 1:nout ) {
    w1 <- rnorm( n1   , 0, sd_w1 )
    b1 <- rnorm( n1   , 0, sd_b1 )
    w2 <- rnorm( n1*n2, 0, sd_w2 )
    b2 <- rnorm(    n2, 0, sd_b2 )
    w3 <- rnorm(    n2, 0, sd_w3 )
    b3 <- rnorm(     1, 0, sd_b3 )
    out[,j] <- fNN.1nn1( x, b1, b2, b3, w1, w2, w3 )
  }
  m <- rowMeans(out) ; s <- apply(out, 1, sd)
  plot  (x, m           , type="l", lwd=3,
                          main=expression("1-10-10-1 (n"[p]*" = 141)"),
                          xlab="", ylab="", ylim=c(-30,30) )
  points(x, m-s         , type="l", lty=2 )
  points(x, m+s         , type="l", lty=2 )
  points(x, out[,1     ], type="l", col="red"   )
  points(x, out[,nout/2], type="l", col="blue"  )
  points(x, out[,nout  ], type="l", col="green" )

  x <- c( -2,  0,    2 )
  y <- c( -10, 10,  -1 ) ; ny <- length(y) ; Sy <- diag(1,ny)
  niNN <- 2e4
  logPrior <- function( b, mb.=mb, Sb.=Sb ) { dmvnorm( b, mb., Sb., log=T ) }
  logLik   <- function( b, y.=y, Sy.=Sy, x.=x, f.=f ) {
    dmvnorm( y., mean=f.(x, b), sigma=Sy., log=T ) }
  logPost  <- function( b, mb.=mb, Sb.=Sb, y.=y, Sy.=Sy, x.=x, f.=f ) {
    logPrior(b, mb., Sb.) + logLik(b, y., Sy., x., f.) }
  MetropolisLogPostMAP <- function( f.=f, logp=logPost, b0=mb,
                                    SProp=Sb/100, ni=1e4, ... ) {
    bChain <- matrix( NA, nrow=ni, ncol=length(b0) ) ; bChain[1,] <- b0
    bMAP   <- b0 ; logpbMAP <- logp(b0,f.=f.,...)
    for( i in 2:ni ) {
       logpb0 <- logp(b0,f.=f.,...)
       b1     <- bChain[i-1,] + as.numeric( rmvnorm(1,sigma=SProp) ) 
       logpb1 <- logp(b1,f.=f.,...)
       if( log(runif(1)) < (logpb1 - logpb0) ) {
         bChain[i,] <- b0 <- b1
         if( logpb1 > logpbMAP ) {bMAP <- b1 ; logpbMAP <- logpb1} }
       else bChain[i,]   <- b0
    }
    return( list( bChain, bMAP ) )
  }

  fNN.xb.1n1 <- function( x, b ) {
    np <- length(b) ; n1 <- ( length(b)-1 ) / 3
    b1 <- b[       1  :  n1      ]
    b2 <- b[  n1  +1             ]
    w1 <- b[ (n1  +2) : (n1*2+1) ]
    w2 <- b[ (n1*2+2) : np       ]
  return( fNN.1n1(x,b1,b2,w1,w2) ) }
  fNN.xb.1nn1 <- function( x, b ) {
    np <- length(b) ; n1 <- n2 <- sqrt( np+3 ) - 2
    b1 <- b[          1  :  n1         ]
    b2 <- b[ (n1     +1) : (n1  +n2  ) ]
    b3 <- b[  n1  +n2+1                ]
    w1 <- b[ (n1  +n2+2) : (n1*2+n2+1) ]
    w2 <- b[ (n1*2+n2+2) : (np  -n2  ) ]
    w3 <- b[ (np  -n2+1) :  np         ]
  return( fNN.1nn1(x,b1,b2,b3,w1,w2,w3) ) }

  mw1 <- mb1 <- mw2 <- mb2 <- 0
  Vw1 <- Vb1 <- Vw2 <- Vb2 <- 9
  set.seed(1)
  n1.A   <- 4
  mbNN.A <-       c( rep(mw1,n1.A), rep(mb1,n1.A), rep(mw2,n1.A), mb2 )
  SbNN.A <- diag( c( rep(Vw1,n1.A), rep(Vb1,n1.A), rep(Vw2,n1.A), Vb2 ) )
  bNNpriorSample.A <- rmvnorm(niNN, mbNN.A, SbNN.A)
  bNNpost.A   <- MetropolisLogPostMAP( fNN.xb.1n1, b0=mbNN.A, SProp=SbNN.A/100,
                                       mb.=mbNN.A, Sb.=SbNN.A, y.=y, ni=niNN*10 )
  bNNChain.A  <- bNNpost.A[[1]] ; fAccepted.A <- dim(unique(bNNChain.A))[1] / niNN
  bNNMAP.A    <- bNNpost.A[[2]] ; bNNmean.A   <- colMeans(bNNChain.A)
  n1.B   <- 20
  mbNN.B <-       c( rep(mw1,n1.B), rep(mb1,n1.B), rep(mw2,n1.B), mb2 )
  SbNN.B <- diag( c( rep(Vw1,n1.B), rep(Vb1,n1.B), rep(Vw2,n1.B), Vb2 ) )
  bNNpriorSample.B <- rmvnorm(niNN, mbNN.B, SbNN.B)
  bNNpost.B  <- MetropolisLogPostMAP( fNN.xb.1n1, b0=mbNN.B, SProp=SbNN.B/100, 
                                      mb.=mbNN.B, Sb.=SbNN.B, y.=y, ni=niNN )
  bNNChain.B <- bNNpost.B[[1]] ; fAccepted.B <- dim(unique(bNNChain.B))[1] / niNN
  bNNMAP.B   <- bNNpost.B[[2]] ; bNNmean.B   <- colMeans(bNNChain.B)

# Marginal posterior distributions for the 13 parameters of a 1-4-1 network.
  n1     <- 4
  pNames <- c(paste0("w1.",1:n1),paste0("b1.",1:n1),paste0("w2.",1:n1),"b2")
  par( mfrow=c(4,4), mar=c(2,2,0,0) )
  for(i in  1: 9){ hist( bNNChain.A[,i], main="", xlab="", ylab="", xlim=c(-12,12), xaxt="n" )
                   # axis(1, labels=FALSE, at=-c(-3,3))
                   legend("topleft",pNames[i]) }
  for(i in 10:13){ hist( bNNChain.A[,i], main="", xlab="", ylab="", xlim=c(-12,12) )
                   legend("topleft",pNames[i]) }

  mw1 <- mb1 <- mw2 <- mb2 <- mw3 <- mb3 <- 0
  Vw1 <- Vb1 <- Vw2 <- Vb2 <- Vw3 <- Vb3 <- 9
  set.seed(1)
  n1.D   <- n2.D <- 2
  mbNN.D <-       c( rep(mb1,n1.D), rep(mb2,     n2.D),     mb3      , 
                     rep(mw1,n1.D), rep(mw2,n1.D*n2.D), rep(mw3,n2.D) )
  SbNN.D <- diag( c( rep(Vb1,n1.D), rep(Vb2,     n2.D),     Vb3      , 
                     rep(Vw1,n1.D), rep(Vw2,n1.D*n2.D), rep(Vw3,n2.D) ) )
  bNNpriorSample.D <- rmvnorm(niNN, mbNN.D, SbNN.D)
  bNNpost.D   <- MetropolisLogPostMAP( fNN.xb.1nn1, b0=mbNN.D, SProp=SbNN.D/100,
                                       mb.=mbNN.D, Sb.=SbNN.D, y.=y, ni=niNN )
  bNNChain.D  <- bNNpost.D[[1]] ; fAccepted.D <- dim(unique(bNNChain.D))[1] / niNN
  bNNMAP.D    <- bNNpost.D[[2]] ; bNNmean.D   <- colMeans(bNNChain.D)
  n1.E   <- n2.E <- 6
  mbNN.E <-       c( rep(mb1,n1.E), rep(mb2,     n2.E),     mb3      , 
                     rep(mw1,n1.E), rep(mw2,n1.E*n2.E), rep(mw3,n2.E) )
  SbNN.E <- diag( c( rep(Vb1,n1.E), rep(Vb2,     n2.E),     Vb3      , 
                     rep(Vw1,n1.E), rep(Vw2,n1.E*n2.E), rep(Vw3,n2.E) ) )
  bNNpriorSample.E <- rmvnorm(niNN, mbNN.E, SbNN.E)
  bNNpost.E   <- MetropolisLogPostMAP( fNN.xb.1nn1, b0=mbNN.E, SProp=SbNN.E/100,
                                       mb.=mbNN.E, Sb.=SbNN.E, y.=y, ni=niNN )
  bNNChain.E  <- bNNpost.E[[1]] ; fAccepted.E <- dim(unique(bNNChain.E))[1] / niNN
  bNNMAP.E    <- bNNpost.E[[2]] ; bNNmean.E   <- colMeans(bNNChain.E)
  n1.F   <- n2.F <- 10
  mbNN.F <-       c( rep(mb1,n1.F), rep(mb2,     n2.F),     mb3      , 
                     rep(mw1,n1.F), rep(mw2,n1.F*n2.F), rep(mw3,n2.F) )
  SbNN.F <- diag( c( rep(Vb1,n1.F), rep(Vb2,     n2.F),     Vb3      , 
                     rep(Vw1,n1.F), rep(Vw2,n1.F*n2.F), rep(Vw3,n2.F) ) )
  bNNpriorSample.F <- rmvnorm(niNN, mbNN.F, SbNN.F)
  bNNpost.F   <- MetropolisLogPostMAP( fNN.xb.1nn1, b0=mbNN.F, SProp=SbNN.F/100,
                                       mb.=mbNN.F, Sb.=SbNN.F, y.=y, ni=niNN )
  bNNChain.F  <- bNNpost.F[[1]] ; fAccepted.F <- dim(unique(bNNChain.F))[1] / niNN
  bNNMAP.F    <- bNNpost.F[[2]] ; bNNmean.F   <- colMeans(bNNChain.F)

# Bayesian calibration of neural networks. Top row: one hidden layer.
# Bottom row: two hidden layers.
# Each line is generated by a different parameter vector: MAP (thick black),
# mean (thick red), random vectors from the MCMC (thin black).
  nxpred <- 500 ; xpred <- seq(-5,5,length.out=nx)
  isubsample <- round( seq(niNN/2,niNN,length.out=10) )
  par(mfrow=c(2,3), mar=c(2,3,2,1))
  plot( xpred, fNN.xb.1n1(xpred,bNNMAP.A), type="l", lwd=3,
        xaxt='n', ylim=c(-20,20), xlab="", ylab="", main="1-4-1" )
  points( xpred, fNN.xb.1n1(xpred,bNNmean.A), type="l", lwd=3, col="red")
  points( cbind(x,y), col='blue', cex=2, lwd=3 )
  for(i in isubsample) { lines( xpred, fNN.xb.1n1(xpred,bNNChain.A[i,]), lwd=0.5 ) }
  plot( xpred, fNN.xb.1n1(xpred,bNNMAP.B), type="l", lwd=3,
        xaxt='n', yaxt='n', ylim=c(-20,20), xlab="x", ylab="", main="1-20-1" )
  points( xpred, fNN.xb.1n1(xpred,bNNmean.B), type="l", lwd=3, col="red")
  points( cbind(x,y), col='blue', cex=2, lwd=3 )
  for(i in isubsample) { lines( xpred, fNN.xb.1n1(xpred,bNNChain.B[i,]),
                                lwd=0.5 ) }
  plot.new()
  legend( "center", legend=c("MAP","mean","sample"),
          col=c("black","red","black"), lty=c(1,1,1), lwd=c(3,3,0.5), cex=1 )
  plot( xpred, fNN.xb.1nn1(xpred,bNNMAP.D), type="l", lwd=3,
        ylim=c(-20,20), xlab="", ylab="", main="1-2-2-1" )
  points( xpred, fNN.xb.1nn1(xpred,bNNmean.D), type="l", lwd=3, col="red")
  points( cbind(x,y), col='blue', cex=2, lwd=3 )
  for(i in isubsample) { lines( xpred, fNN.xb.1nn1(xpred,bNNChain.D[i,]), lwd=0.5 ) }
  plot( xpred, fNN.xb.1nn1(xpred,bNNMAP.E), type="l", lwd=3,
        yaxt='n', ylim=c(-20,20), xlab="x", ylab="", main="1-6-6-1" )
  points( xpred, fNN.xb.1nn1(xpred,bNNmean.E), type="l", lwd=3, col="red")
  points( cbind(x,y), col='blue', cex=2, lwd=3 )
  for(i in isubsample) { lines( xpred, fNN.xb.1nn1(xpred,bNNChain.E[i,]),
                                lwd=0.5 ) }
  plot( xpred, fNN.xb.1nn1(xpred,bNNMAP.F), type="l", lwd=3,
        yaxt='n', ylim=c(-20,20), xlab="x", ylab="", main="1-10-10-1" )
  points( xpred, fNN.xb.1nn1(xpred,bNNmean.F), type="l", lwd=3, col="red")
  points( cbind(x,y), col='blue', cex=2, lwd=3 )
  for(i in isubsample) { lines( xpred, fNN.xb.1nn1(xpred,bNNChain.F[i,]),
                                lwd=0.5 ) }

  k <- 0.1 ; SbNN.B2 <- SbNN.B * k
  set.seed(1)
  bNNpriorSample.B2 <- rmvnorm(niNN, mbNN.B, SbNN.B2)
  bNNpost.B2  <- MetropolisLogPostMAP( fNN.xb.1n1, b0=mbNN.B, SProp=SbNN.B2/100, 
                                     mb.=mbNN.B, Sb.=SbNN.B2, y.=y, ni=niNN )
  bNNChain.B2 <- bNNpost.B2[[1]] ; fAccepted.B2 <- dim(unique(bNNChain.B2))[1] / niNN
  bNNMAP.B2   <- bNNpost.B2[[2]] ; bNNmean.B2   <- colMeans  (bNNChain.B2)

# Prior and posterior predictions of a 20-neuron one-hidden-layer neural network
# with narrow prior for weights and biases.
  nxpred <- 500 ; xpred <- seq(-5,5,length.out=nx)
  isubsample <- round( seq(niNN/2,niNN,length.out=10) )
  par( mfrow=c(1,2), mar=c(2,3,2,1) )
  plot( xpred, fNN.xb.1n1(xpred,mbNN.B), type="l", lwd=3,
        ylim=c(-20,20), xlab="", ylab="", main="1-20-1 (prior)", col="red" )
  for(i in isubsample) { lines( xpred, fNN.xb.1n1(xpred,bNNpriorSample.B2[i,]),
                                lwd=0.5 ) }
  legend( "topleft",
          legend=c("mean","sample"),
          col=c("red","black"), lty=c(1,1), lwd=c(3,0.5), cex=0.5 )
  plot( xpred, fNN.xb.1n1(xpred,bNNMAP.B2), type="l", lwd=3, yaxt="n",
        ylim=c(-20,20), xlab="", ylab="", main="1-20-1 (posterior)" )
  points( xpred, fNN.xb.1n1(xpred,bNNmean.B2), type="l", lwd=3, col="red")
  points( cbind(x,y), col='blue', cex=2, lwd=3 )
  for(i in isubsample) { lines( xpred, fNN.xb.1n1(xpred,bNNChain.B2[i,]),
                                lwd=0.5 ) }
  legend( "topleft",
          legend=c("MAP","mean","sample"),
          col=c("black","red","black"), lty=c(1,1,1), lwd=c(3,3,0.5), cex=0.5 )

  mnist <- dslabs::read_mnist( download=T, destdir=getwd() )
  # mnist <- dslabs::read_mnist( getwd() )

  blockSqMat = function( mat.in, block=2 ) {
    size.in <- nrow(mat.in) ; size.out <- size.in / block
    mat.cmn <- matrix( rep(NA,size.in*size.out), nrow=size.in  )
    mat.out <- matrix( rep(NA,size.out^2      ), nrow=size.out )
    for(c in 1:size.out) {
      c.in        <- seq( (c-1)*block+1, c*block )
      mat.cmn[,c] <- apply( mat.in[,c.in], 1, mean )    
    } # Average columns by block
    for(r in 1:size.out) {
      r.in        <- seq( (r-1)*block+1, r*block )
      mat.out[r,] <- apply( mat.cmn[r.in,], 2, mean )    
    } # Average rows by block
    return(mat.out)
  }

  block <- 14 ; size.out <- 28 / block
  nim.t <- 100 ; mnist$train$images.bl <- matrix( NA, nim.t, size.out^2 )
  for(i in 1:nim.t) {  mtrain <- matrix( mnist$train$images[i,], nrow=28 )
                    mtrain.bl <- blockSqMat( mtrain, block )
                  vectrain.bl <- mtrain.bl[ 1 : size.out^2 ]
    mnist$train$images.bl[i,] <- vectrain.bl }
  nim.u <- 100 ; mnist$test$images.bl  <- matrix( NA, nim.u, size.out^2 )
  for(i in 1:nim.u) {   mtest <- matrix( mnist$test$images [i,], nrow=28 )
                     mtest.bl <- blockSqMat( mtest , block )
                   vectest.bl <- mtest.bl [ 1 : size.out^2 ]
     mnist$test$images.bl[i,] <- vectest.bl }

  stsoftmax <- function(y){ exp(y) / sum(exp(y)) }

  fNN.4310 <- function( x, b1, b2, w1, w2, fh=tanh, fy=stsoftmax ) {
    # nx <- length(x) ; nh <- length(b1) ; ny <- length(b2)
    nh     <- 3 ; ny <- 10
    w1.1   <- w1[      1 : nh   ] ; w1.2 <- w1[(nh  +1):(nh*2)]
    w1.3   <- w1[(nh*2+1):(nh*3)] ; w1.4 <- w1[(nh*3+1):(nh*4)]
    w2.1   <- w2[      1 : ny   ] ; w2.2 <- w2[(ny  +1):(ny*2)]
    w2.3   <- w2[(ny*2+1):(ny*3)]
    X      <- c    ( 1, x )                           # [1 x  5]
    theta1 <- rbind( b1, w1.1, w1.2, w1.3, w1.4 )     # [5 x  3]
    theta2 <- rbind( b2, w2.1, w2.2, w2.3 )           # [4 x 10]
    h      <- fh( X %*% theta1 ) ; H <- cbind( 1, h ) # [1 x  3] ; [1 x 4]
    y      <- fy( H %*% theta2 )                      # [1 x 10]
    return( y )
  }

  fNN.xb.4310 <- function( x, b ) {
    b1 <- b[1:3] ; b2 <- b[4:13] ; w1 <- b[14:25] ; w2 <- b[26:55]
    nx <- length(x) / 4 ; if(nx==1) x <- t(x) # Ensure number of x-columns = 4
    y  <- matrix(NA,nrow=nx,ncol=10)
    for(r in 1:nx) y[r,] <- fNN.4310( x[r,], b1, b2, w1, w2 )
    return( y )
  }

  x.t  <- mnist$train$images.bl ; xmax <- max(x.t)  
  x.ts <- x.t / xmax # Scaling the inputs to [0,1]
  x.u  <- mnist$test$images.bl
  x.us <- x.u / xmax # Scaling the inputs as we did for the training data

  labels.t <- mnist$train$labels[1:nim.t]
  y.t      <- matrix( NA, nrow=nim.t, ncol=10 )
  for(r in which(labels.t==0)) { y.t[r,] <- c(1,0,0,0,0,0,0,0,0,0) }
  for(r in which(labels.t==1)) { y.t[r,] <- c(0,1,0,0,0,0,0,0,0,0) }
  for(r in which(labels.t==2)) { y.t[r,] <- c(0,0,1,0,0,0,0,0,0,0) }
  for(r in which(labels.t==3)) { y.t[r,] <- c(0,0,0,1,0,0,0,0,0,0) }
  for(r in which(labels.t==4)) { y.t[r,] <- c(0,0,0,0,1,0,0,0,0,0) }
  for(r in which(labels.t==5)) { y.t[r,] <- c(0,0,0,0,0,1,0,0,0,0) }
  for(r in which(labels.t==6)) { y.t[r,] <- c(0,0,0,0,0,0,1,0,0,0) }
  for(r in which(labels.t==7)) { y.t[r,] <- c(0,0,0,0,0,0,0,1,0,0) }
  for(r in which(labels.t==8)) { y.t[r,] <- c(0,0,0,0,0,0,0,0,1,0) }
  for(r in which(labels.t==9)) { y.t[r,] <- c(0,0,0,0,0,0,0,0,0,1) }
  labels.u <- mnist$test$labels[1:nim.u]
  y.u      <- matrix(NA,nrow=nim.u,ncol=10)
  for(r in which(labels.u==0)) { y.u[r,] <- c(1,0,0,0,0,0,0,0,0,0) }
  for(r in which(labels.u==1)) { y.u[r,] <- c(0,1,0,0,0,0,0,0,0,0) }
  for(r in which(labels.u==2)) { y.u[r,] <- c(0,0,1,0,0,0,0,0,0,0) }
  for(r in which(labels.u==3)) { y.u[r,] <- c(0,0,0,1,0,0,0,0,0,0) }
  for(r in which(labels.u==4)) { y.u[r,] <- c(0,0,0,0,1,0,0,0,0,0) }
  for(r in which(labels.u==5)) { y.u[r,] <- c(0,0,0,0,0,1,0,0,0,0) }
  for(r in which(labels.u==6)) { y.u[r,] <- c(0,0,0,0,0,0,1,0,0,0) }
  for(r in which(labels.u==7)) { y.u[r,] <- c(0,0,0,0,0,0,0,1,0,0) }
  for(r in which(labels.u==8)) { y.u[r,] <- c(0,0,0,0,0,0,0,0,1,0) }
  for(r in which(labels.u==9)) { y.u[r,] <- c(0,0,0,0,0,0,0,0,0,1) }

# First four MNIST training matrices. Top row: original data.
# Bottom row: block-averaged. From left to right, the digits are 5, 0, 4, 1.
  layout( matrix(1:8, 2, 4, byrow=T) )
  par( mar=c(1,1,1,1) , omi=c(0,0,0.5,0) )
  for(i in 1:4) {
    mati <- matrix( mnist$train$images[i,], nrow=28 )
    labl <- mnist$train$labels[i]
    image(1:28, 1:28, mati[,28:1], main="", xaxt="n", yaxt="n", bty="n",
      zlim=c(0,255), col=gray(seq(1,0,-0.05)), xlab="", ylab="", asp=1 )
  }
  for(i in 1:4) {
    mati <- matrix( x.t[i,], nrow=size.out ) ; labl <- labels.t[i]
    image(1:size.out, 1:size.out, mati[,size.out:1], main="", xaxt="n", yaxt="n",
          zlim=c(0,255), bty="n", col=gray(seq(1,0,-0.05)), xlab="", ylab="", asp=1 )
  }

  mbNN4310 <-       c( rep(0,3), rep(0,10), rep(0.5,12), rep(1,30) )
  SbNN4310 <- diag( c( rep(1,3), rep(1,10), rep(1  ,12), rep(9,30) ) )
  logLik.Cl <- function( b, y.=y, x.=x, f.=f ) {
    sum( y.*log(f.(x.,b)) + (1-y.)*log(1-f.(x.,b)) ) }
  logPost.Cl <- function( b, mb.=mb, Sb.=Sb, y.=y, x.=x, f.=f ) {
    logPrior(b, mb., Sb.) + logLik.Cl(b, y., x., f.) }
  set.seed(1) ; ni.t <- 1e4
  bNN4310post  <- MetropolisLogPostMAP( f.=fNN.xb.4310, logp=logPost.Cl,
    b0=mbNN4310, SProp=SbNN4310/1000, mb.=mbNN4310, Sb.=SbNN4310,
    y.=y.t, x.=x.ts, ni=ni.t )
  bNN4310Chain <- bNN4310post[[1]]
  fAccepted    <- dim(unique(bNN4310Chain))[1] / ni.t ; fAccepted
  bNN4310MAP   <- bNN4310post[[2]] ; bNN4310mean <- colMeans(bNN4310Chain)
  theta <- bNN4310MAP 
  c( logPrior  ( theta, mbNN4310, SbNN4310 ),
     logLik.Cl ( theta, y.t, x.ts, fNN.xb.4310 ),
     logPost.Cl( theta, mbNN4310, SbNN4310, y.t, x.ts, fNN.xb.4310 ) )

  yMAP.t <- fNN.xb.4310( x.ts, bNN4310MAP ) ; round( yMAP.t, 2 )
  nt       <- length(labels.t) # = dim(y.t)[1]
  classMAP.t <- rep( NA, nt )
  for(r in 1:nt) { classMAP.t[r] <- which.max(yMAP.t[r,]) - 1 }
  fCorrect <- sum( classMAP.t==labels.t ) / nt ; fCorrect

  yMAP.u     <- fNN.xb.4310( x.us, bNN4310MAP ) ; round( yMAP.u, 2 )
  nu         <- length(labels.u) # = dim(y.t)[1]
  classMAP.u <- rep( NA, nu )
  for(r in 1:nu) { classMAP.u[r] <- which.max(yMAP.u[r,]) - 1 }
  fCorrect   <- sum( classMAP.u==labels.u ) / nu ; fCorrect

  n.smpl <- 20
  isubsample <- round( seq(ni.t/2,ni.t,length.out=n.smpl) )
  fCorrect.sample <- rep( NA, n.smpl )
  for(i in 1:n.smpl) {
    bNN4310.i  <- bNN4310Chain[ isubsample[i], ]
    y.u.i      <- fNN.xb.4310( x.us, bNN4310.i )
    class.u.i  <- rep( NA, nu )
    for(r in 1:nu) { class.u.i[r] <- which.max(y.u.i[r,]) - 1 }
    fCorrect.sample[i] <- sum( class.u.i==labels.u ) / nu
  }
  range(fCorrect.sample)

# Classification of MNIST digit-images using a 4-3-10 network with the MAP
# parameter vector from Bayesian calibration. Panel headers indicate true
# digit-values. Bars represent the network's average probability of assigning
# the digits 0..9, in each case distinguishing correct (red) and incorrect
# attributions (black).
  layout( matrix(1:10, 2, 5, byrow=T) )
  par( mar=c(3,3,3,1), omi=c(0,0,0,0) )
  for(i in 0) {
    pMAP.i <- colMeans( yMAP.u[labels.u==i,] )
    cols   <- rep("black",10) ; cols[i+1] <- "red"
    barplot( pMAP.i , main=i, col=cols, ylim=c(0,0.45) ) }
  for(i in 1:4) {
    pMAP.i <- colMeans( yMAP.u[labels.u==i,] )
    cols   <- rep("black",10) ; cols[i+1] <- "red"
    barplot( pMAP.i , main=i, col=cols, ylim=c(0,0.45), yaxt="n" ) }
  for(i in 5) {
    pMAP.i <- colMeans( yMAP.u[labels.u==i,] )
    cols   <- rep("black",10) ; cols[i+1] <- "red"
    barplot( pMAP.i , main=i, col=cols, ylim=c(0,0.45), names.arg=0:9 ) }
  for(i in 6:9) {
    pMAP.i <- colMeans( yMAP.u[labels.u==i,] )
    cols   <- rep("black",10) ; cols[i+1] <- "red"
    barplot( pMAP.i , main=i, col=cols, ylim=c(0,0.45), yaxt="n", names.arg=0:9 ) }

  n <- 1e4
  pi.url <- "http://oeis.org/A000796/b000796.txt"
  pi.dec <- read.csv( pi.url, sep=" ")$X3[ 1:n ]
  n.t        <- n/2
  pi.dec.t   <- pi.dec  [1:n.t] ; pi.dec.u <- pi.dec  [(n.t+1):n]
  matPtrans <- function( seq ) {
    lab   <- sort( unique(seq) ) ; n.lab <- length( lab )
    f.mat <- matrix( NA, n.lab, n.lab, dimnames=list(lab,lab) )
    for(r in 1:n.lab) {
      for(c in 1:n.lab){ f.mat[r,c] <- length( which(
        seq[1:(n-1)]==lab[r] & seq[2: n   ]==lab[c] ) )
      }
    }
    f.int <- rowSums(f.mat)
    p.mat <- f.mat / f.int ; round( p.mat, 2 )
    return( list( L=lab, P=round(p.mat,2) ) )
  }

  print( matPtrans( pi.dec.t ) [["P"]] [1:3,1:10] )
  cat( "..." )

# Bertrand Russell - The Problems Of Philosophy
  f         <- "https://www.gutenberg.org/cache/epub/5827/pg5827.txt"
  txt       <- scan( f, "char", skip=66, blank.lines.skip=T, nlines=1000 )
  txt       <- gsub( "\\?"        , "\\.", txt )
  txt       <- gsub( "\\_"        , ""   , txt )
  txt       <- gsub( "([,;:()])"  , ""   , txt )
  txt       <- gsub( "--i.e."     , ""   , txt )
  txt       <- gsub("[[:digit:]]+", ""   , txt )
  txt       <- gsub( "-"          , " "  , txt )
  txt       <- gsub( "\n"         , " "  , txt )
  txt       <- gsub( "  "         , " "  , txt )
  txt       <- txt[ which(txt!="") ]
    ifin    <- grep("\\w+\\.", txt) ; nsen <- length(ifin)
  txt[ifin] <- substr( txt[ifin], 1, nchar(txt[ifin])-1 )
    itxt    <- c( seq_along(txt), ifin+0.5 )
  txt       <- c( txt, rep(".",nsen) )
  txt       <- txt[ order(itxt) ] ; n.txt <- length(txt)
  txt       <- unlist( lapply( as.list(txt),
                               function(x) strsplit(x,split=" ") ) )
  fwords <- table(txt)           ; n.words <- length(fwords)
  pwords <- fwords / sum(fwords) ;   words <- names (fwords)

# Vocabulary of the first 1000 lines of Russell (2012).
# The 1552 different words are sorted alphabetically from left to right.
# Vertical position is log10( word frequency ).
  par(mfrow=c(1,1))
  txt.num          <- match( txt, words )
  hist.data        <- hist( txt.num, breaks=0.5:(n.words+0.5), plot=FALSE )
  hist.data$counts <- log10(hist.data$counts)
  plot( 1:n.words, hist.data$counts, ylab='log10(frequency)',
        frame=FALSE, pch=".", xaxt="n", xlab="", cex.axis=0.5, cex.lab=0.5 )
  text( x=1:n.words, y=hist.data$counts, labels=words, cex=0.5 )  

  nextnum <- function( num, v=txt.num ) {
    nextnum <- v[ 1 + which(v==num) ]
    nextnum <- nextnum[ !is.na(nextnum) ]
    return( nextnum ) }

  rw1 <- function( num, v=txt.num ) { 
    sample( nextnum(num,v), 1, replace=TRUE ) }

  set.seed(7)
  prompt      <- "." ; lastnum <- match( prompt, words )
  txt.num.out <- NULL
  for(i in 1:100) {
    txt.num.out <- c( txt.num.out, rw1(lastnum) )
    lastnum     <- tail( txt.num.out, 1 ) }
  txt.out <- words[ txt.num.out ]

  n       <- 150 ; nz <- 3 ; ny <- 4
  z       <- as.numeric(iris[,5])  ; y      <- as.matrix(iris[,1:ny])
  r.train <- c(1:25,51:75,101:125) ; r.test <- setdiff(1:150,r.train)
  y.train <- y[ r.train, ]         ; y.test <- y[ r.test, ]
  z.train <- z[ r.train  ]         ; z.test <- z[ r.test  ]
  m.train <- s.train <- vector("list",nz)
  for( j in 1:nz ) {
    m.train[[j]] <- apply( y.train[z.train==j,], 2, mean )
    s.train[[j]] <- apply( y.train[z.train==j,], 2, sd   ) }

  L.iris <- function( y, z, m=m.train, s=s.train ) {
    prod( dnorm( y, m[[z]], s[[z]] ) ) }
  n.test <- length(z.test)
  L.test <- pz_y.test <- matrix( NA, nrow=n.test, ncol=nz )
  for(j in 1:nz) {
    L.test[,j] <- apply( y.test, 1, function(y) { L.iris( y, z=j ) } ) }
  pz_y.test  <- t( apply( L.test, 1, function(Ly) { Ly / sum(Ly) } ) )
  class.test <- apply( L.test, 1, which.max )
  accuracy   <- sum( class.test==z.test) / n.test ; accuracy
  pclass_y.test <- pz_y.test  [ cbind(1:n.test,class.test  ) ]
  mean( pclass_y.test )
  ierror        <- which( class.test != z.test )
  cbind( ierror, z.test[ierror  ], pz_y.test[ierror,] )

  S.train <- vector("list",nz)
  for( j in 1:nz ) { S.train[[j]] <- cov( y.train[z.train==j,] ) }
  L.iris.S <- function( y, z, m=m.train, S=S.train ) {
    dmvnorm( y, m[[z]], S[[z]] ) }
  L.test.S <- pz_y.test.S <- matrix( NA, nrow=n.test, ncol=nz )
  for(j in 1:nz) {
    L.test.S[,j] <- apply( y.test, 1, function(y) { L.iris.S( y, z=j ) } ) }
  pz_y.test.S  <- t( apply( L.test.S, 1, function(Ly) { Ly / sum(Ly) } ) )
  class.test.S <- apply( L.test.S, 1, which.max )
  accuracy.S   <- sum( class.test.S==z.test) / n.test ; accuracy.S
  pclass_y.test.S <- pz_y.test.S[ cbind(1:n.test,class.test.S) ]
  mean( pclass_y.test.S )
  ierror.S <- which( class.test.S != z.test )
  cbind( ierror.S, z.test[ierror.S], pz_y.test.S[ierror.S,] )

# Gaussian classification of iris flowers.
# Left: covariances ignored. Right: covariances accounted for.
# Digits 1 to 3 indicate the most probable species for each of 75 tested plants.
# Classification is correct (black digits) or incorrect (large red digits).
  par( mfrow=c(1,2), mar=c(4,4,3,1) )
  cols   <- rep("black",n.test) ; cols  [ierror  ] <- "red"
  cols.S <- rep("black",n.test) ; cols.S[ierror.S] <- "red"
  cexs   <- rep( 0.5,   n.test) ; cexs  [ierror  ] <- 0.75
  cexs.S <- rep( 0.5,   n.test) ; cexs.S[ierror.S] <- 0.75
  plot( pclass_y.test, pch=as.character(class.test), col=cols,
        ylim=range(pclass_y.test, pclass_y.test.S),
        main="Naive Bayes", cex=cexs, ylab="max( p[z] )")
        # main="Naive Bayes", cex=0.5, ylab="max( p[z] )")
  plot( pclass_y.test.S, pch=as.character(class.test.S), col=cols.S,
        ylim=range(pclass_y.test, pclass_y.test.S), yaxt="n", ylab="",
        main="Bayes", cex=cexs.S)
