## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 14. Thirteen Ways to Fit a Straight Line

  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(mvtnorm)
  library(rjags)
  library(rsvg)

  x <- c(10,20,30) ; y <- c(6.09,8.81,10.66) ; ny <- length(y)
  
  b2 <- cov(x,y) / var(x)
  b1 <- mean(y) - b2 * mean(x)

  lm.y <- lm(y~x) ; print(coef(lm.y))
# Linear regression on three data points.
  plot(x,y) ; abline(lm.y)

  Vy   <- rep(3,ny) ; Sy <- diag(Vy)
  X    <- cbind( 1, x )                       # Design matrix
  Sb_y <- solve(   t(X) %*% solve(Sy) %*% X ) # GLS variance
  mb_y <- Sb_y %*% t(X) %*% solve(Sy) %*% y   # GLS estimator

  mb <- c(0   ,0   ) ; nb <- length(mb)
  Vb <- c(1.e4,1.e4) ; Sb <- diag(Vb)
  
  Sb_y_LS72 <- solve( solve(Sb) + t(X) %*% solve(Sy) %*% X ) ; 
  mb_y_LS72 <- Sb_y_LS72 %*% (solve(Sb) %*% mb + t(X) %*% solve(Sy) %*% y)

  KalmanGain <- function(B,H,R){ B %*% t(H) %*% solve(H %*% B %*% t(H) + R) }
  K          <- KalmanGain( Sb, X, Sy )
  mb_y       <- mb + K %*% ( y - X %*% mb )
  Sb_y       <- ( diag(nb) - K %*% X ) %*% Sb

  GaussCond <- function( mz, Sz, y ) {
    i <- 1 : ( length(mz) - length(y) )
    m  <- mz[i]   + Sz[i,-i] %*% solve(Sz[-i,-i]) %*% (y-mz[-i])
    S  <- Sz[i,i] - Sz[i,-i] %*% solve(Sz[-i,-i]) %*% Sz[-i,i]
    return( list( m=m, S=S ) ) }
  mz <- c( mb, X %*% mb )
  Sz <- rbind( cbind( Sb        ,   Vb*t(X)            ),
               cbind( t(Vb*t(X)), t(Vb*t(X))%*%t(X)+Sy ) )
  Post <- GaussCond( mz, Sz, y ) ; mb_y <- Post$m ; Sb_y <- Post$S

  DAGlm <- grViz( "digraph{ graph[ rankdir=BT, ranksep=.5, nodesep=.3 ]
    node[ shape=rectangle ]
      b1[label='@@1'] ; b2[label='@@2']
    subgraph cluster_0 {
      graph[ shape=rectangle ]
        style=rounded ; fontsize=20 ; labelloc=b ; label='Data'
      node[ shape = ellipse ]
        y1[label='y[1]'] ; y2[label='y[2]'] ; y3[label='y[3]'] }
    edge[ arrowhead=vee, arrowsize=1, color=black ]
      b1 -> {y1,y2,y3} [label='1']
    edge[ arrowhead=vee, arrowsize=1, color=red, fontcolor=red ]
      b2 -> y1[label='x[1]'] ; b2 -> y2[label='x[2]'] ; b2 -> y3[label='x[3]'] }
    [1]: paste0('&beta;&#x2081;','\\n(intercept)')
    [2]: paste0('&beta;&#x2082;','\\n(slope)')
  ")
  export_svg(DAGlm) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGlm.png")


# A minimal neural network with a single hidden layer.
  DAGNNLR <- grViz("digraph{ graph[rankdir=LR, ranksep=1, nodesep=0.5]
      subgraph cluster_0 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='OUTPUT']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        y[label = '@@1']
      }
      subgraph cluster_1 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='HIDDEN LAYER']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        h[label = '@@2']
      }
      subgraph cluster_2 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='INPUT']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        x[label = '@@3']
      }
      edge[arrowhead=vee, arrowsize=1.25]
      x -> h -> y }
      [1]: paste0('y')
      [2]: paste0('h')
      [3]: paste0('x')
    ")
  export_svg(DAGNNLR) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGNNLR.png")


  GaussMult <- function( m=rep(0,2), V=rep(1,2) ) {
    n       <- length(m)
    weights <- sapply( 1:n, function(i){ prod(V[-i]) } )
    m_Mult  <- sum( m * weights ) / sum( weights )
    V_Mult  <- prod( V ) / sum( weights )
    return( list( m=m_Mult, V=V_Mult ) ) }
  i  <- 1 ; sy=sqrt(Vy) ; b1 <- 0 ; b2 <- 0
  L <- dnorm( b1 + b2*x[i] - y[i], mean= 0              , sd=sy[i]      )
  L <- dnorm(      b2            , mean=(y[i] - b1)/x[i], sd=sy[i]/x[i] ) / x[i]
  L <- dnorm( b1                 , mean= y[i] - b2 *x[i], sd=sy[i]      )
  ni <- 1e3 ; bChain <- matrix(NA, nrow=ni, ncol=2)
  b1 <- mb[1] ; b2 <- mb[2] ;  bChain[1,] <- c( b1, b2 )
  for (i in (2 : ni)) {
    mb1_L <- y - b2*x ; Vb1_L <- diag( Sy )
    Post1 <- GaussMult( m=c(mb[1],mb1_L), V=c(Vb[1],Vb1_L) )
    mb1   <- Post1$m ; Vb1 <- Post1$V
    b1    <- rnorm( 1, mean=mb1, sd=sqrt(Vb1) )
    mb2_L <- ( y - b1 ) / x ; Vb2_L <- diag( Sy/(x^2) )
    Post2 <- GaussMult( m=c(mb[2],mb2_L), V=c(Vb[2],Vb2_L) )
    mb2   <- Post2$m ; Vb2 <- Post2$V
    b2    <- rnorm( 1, mean=mb2, sd=sqrt(Vb2) )
    bChain[i,] <- c( b1, b2 ) }
  
# Gibbs sampling using own R-code.
  par( mfrow=c(2,2), mar=c(2,2,3,2) )
  plot(bChain[,1],xlab="",ylab="",main="Intercept")
  plot(bChain[,2],xlab="",ylab="",main="Slope")
  hist(bChain[,1],xlab="",ylab="",main="Intercept")
  hist(bChain[,2],xlab="",ylab="",main="Slope")
  
  f.JAGS <- " model { for (i in 1:n) {
      y[i]      ~ dnorm( y.hat[i], tau.y[i] )
      y.hat[i] <- b0 + b1 * x[i]
      tau.y[i] <- pow( sy[i], -1 ) }
    b0 ~ dnorm( 0, .0001 ) ; b1 ~ dnorm( 0, .0001 ) } "
  writeLines( f.JAGS, con="f.JAGS.txt" )
  f.JAGS.data   <- list ( n=ny, y=y, x=x, sy=diag(Sy) )
  f.JAGS.params <- c( "b0", "b1" )
  f.JAGS        <- jags.model( "f.JAGS.txt", data=f.JAGS.data, n.adapt=5e3 )
  update( f.JAGS, n.iter=5e3 )
  f.JAGS.codaSamples <- coda.samples( f.JAGS, var=f.JAGS.params, n.iter=5e3 )
  f.JAGS.mcmcChain   <- as.matrix( f.JAGS.codaSamples )
# Gibbs sampling using JAGS."}
  par( mar=c(2,2,2,2))
  plot( f.JAGS.codaSamples )

# Exercise 3. Kalman Filtering (KF) without intercept uncertainty.
  # K    <- KalmanGain( Sb*diag(c(0,1)), X, Sy )
  # Sb_y <- ( diag(nb) - K %*% X ) %*% Sb*diag(c(0,1))
  # mb_y <- mb + K %*% ( y - X %*% mb )
  # print(round(mb_y,3)) ; print(round(Sb_y,3))