## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 16. Graphical Modelling

  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(mvtnorm)
  library(rsvg)

# Simplified scheme of graphical modelling categories with only the major types
# shown.
  DAGGM <- grViz( "digraph{ graph[]
    node[ shape=box ]
      GM [label='@@1'] ; MRF[label='@@2'] ; MIX[label='@@3']
      BN[label='@@4'] ; GBN[label='@@5'] ; DISCR[label='@@6']
    edge[]
      GM -> {MRF, MIX, BN}
      BN -> {GBN, DISCR} }
    [1]: paste0( 'Graphical models', '\\n(GM)' )
    [2]: paste0( 'Undirected graphs', '\\n(MRF)' )
    [3]: paste0( 'Mixed graphs', '\\n ' )
    [4]: paste0( 'Directed graphs', '\\n(BN)' )
    [5]: paste0( 'Continuous networks', '\\n(GBN)' )
    [6]: paste0( 'Discrete networks', '\\n ' )
  ")
  export_svg(DAGGM) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGGM.png")

# The simplest Directed Acyclic Graph (DAG).
  DAGab <- grViz( "digraph{ graph[rankdir=LR] edge[] A->B }")
  export_svg(DAGab) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGab.png")

# C is conditionally independent of A.
  DAGabc <- grViz( "digraph{ graph[rankdir=LR] edge[] A->B->C }")
  export_svg(DAGabc) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGabc.png")

# DAG with means and conditional variances specified in the node-ellipses and
# edge labelled with the regression coefficient.
  DAGabLabelled <- grViz( "digraph{ graph[rankdir=LR]
    node[]
      A[label=<A<br/>m<SUB>A</SUB>= 0; V<SUB>A</SUB>= 1>]
      B[label=<B<br/>m<SUB>B</SUB>= 0; V<SUB>B</SUB>= 0.75>]
    edge[] A->B [label='r@_{AB}= 0.5'] }
  ")
  export_svg(DAGabLabelled) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGabLabelled.png")

  precMatrix <- function( nodes, Vcond, R ){
    n <- length(nodes) ; W <- 1 / Vcond[1]
    for(k in 2:n){
      rk       <- R[ k, 1:(k-1) ]
      W_top    <- cbind( W * Vcond[k] + rk %*% t(rk), -rk ) / Vcond[k] 
      W_bottom <- cbind(                      -t(rk),   1 ) / Vcond[k]
      W        <- rbind( W_top, W_bottom ) }
    rownames(W) <- colnames(W) <- nodes
    return(W) }

  nodes <- c("A","B") ; Vcond <- c(1,0.75) ; R <- matrix( c(0,0.5,0,0), nrow=2)
  W     <- precMatrix( nodes, Vcond, R )
  Sigma <- solve(W)

  VcondR <- function( W ) {
    n  <- dim(W)[1] ; Vcond <- rep( NA, n )
    Wk <- W         ; R     <- matrix( 0, nrow=n, ncol=n )
    for(k in n:2){
      ik       <- 1 : (k-1)
      Vcond[k] <- 1 / Wk[ k, k]
      R[k,ik]  <-    -Wk[ k,ik] * Vcond[k]
      Wk       <-     Wk[ik,ik] - as.matrix(R[k,ik]) %*% R[k,ik] / Vcond[k] }
    Vcond[1] <- 1 / Wk[1,1]
    return( list( Vcond=Vcond, R=zapsmall(R) ) ) }

  VR <- VcondR(W)

  nodes  <- c("Rain","Irrigation","SoilWater")
  m      <- c( 2    , 0.5        , 1         )
  Vcond  <- c( 0.2  , 0.05       , 0.1       )
  R      <- matrix( rep(0,9), nrow=3) ; R[3,1] <- 0.9 ; R[3,2] <- 1
  W      <- precMatrix(nodes,Vcond,R) ; S <- zapsmall( solve(W) )
  i2     <- 3:1 ; nodes2 <- nodes[i2] ; m2 <- m[i2] ; S2 <- S[i2,i2]
  W2     <- solve(S2) ; VcondR2 <- VcondR(W2)
  Vcond2 <- signif(VcondR2$Vcond,2) ; R2 <- signif(VcondR2$R,2)

# Two DAGS for the same probability distribution. The causal DAG on the left is
# smaller than the non-causal one on the right.
  DAGris <- grViz( "digraph{ graph[]
    node[]
      Rain [label='@@1'] ; Irri [label='@@2'] ; Soil [label='@@3']
      Rain2[label='@@5'] ; Irri2[label='@@6'] ; Soil2[label='@@7']
    edge[]
      Rain  -> Soil [label='@@4-1'] ; Irri  -> Soil  [label='@@4-2']
      Soil2 -> Irri2[label='@@8-1'] ; Soil2 -> Rain2 [label='@@8-2']
      Irri2 -> Rain2[label='@@8-3'] }
    [1]: paste0('Rain'      ,'\\nm=',m[1],',V=',Vcond[1])
    [2]: paste0('Irrigation','\\nm=',m[2],',V=',Vcond[2])
    [3]: paste0('Soil Water','\\nm=',m[3],',V=',Vcond[3])
    [4]: c( paste0('r=',R[3,1]), paste0('r=',R[3,2]) )
    [5]: paste0('Rain'      ,'\\nm=',m2[3],',V=',Vcond2[3])
    [6]: paste0('Irrigation','\\nm=',m2[2],',V=',Vcond2[2])
    [7]: paste0('Soil Water','\\nm=',m2[1],',V=',Vcond2[1])
    [8]: c( paste0('r=',R2[2,1]), paste0('r=',R2[3,1]), paste0('r=',R2[3,2]) )
  ")
  export_svg(DAGris) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGris.png")

  GaussCond <- function( mz, Sz, y ) {
    i <- 1 : ( length(mz) - length(y) )
    m  <- mz[i]   + Sz[i,-i] %*% solve(Sz[-i,-i]) %*% (y-mz[-i])
    S  <- Sz[i,i] - Sz[i,-i] %*% solve(Sz[-i,-i]) %*% Sz[-i,i]
    return( list( m=m, S=S ) ) }

  nodes <- c( "A", "B", "C", "D" )
  m     <- c(  3 ,  4 ,  9 ,  14 )
  S     <- matrix( c(4 ,  4,  8, 12,
                     4 ,  5,  8, 13,
                     8 ,  8, 20, 28,
                     12, 13, 28, 42), nrow=4 )
  VR <- VcondR( solve(S) ) ; Vcond <- signif(VR$Vcond,2) ; R <- VR$R

# DAG with four nodes.
  DAG4 <- grViz( "digraph{ graph[]
    node[]
      A[label='@@1'] ; B[label='@@2'] ; C[label='@@3'] ; D[label='@@4']
    edge[]
      A -> B [label='@@5'] ; A -> C [label='@@6']
      B -> D [label='@@7'] ; C -> D [label='@@8'] }
    [1]: paste0('A','\\nm=',m[1],',V=',Vcond[1])
    [2]: paste0('B','\\nm=',m[2],',V=',Vcond[2])
    [3]: paste0('C','\\nm=',m[3],',V=',Vcond[3])
    [4]: paste0('D','\\nm=',m[4],',V=',Vcond[4])
    [5]: paste0('r=',R[2,1])
    [6]: paste0('r=',R[3,1])
    [7]: paste0('r=',R[4,2])
    [8]: paste0('r=',R[4,3])
  ")
  export_svg(DAG4) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAG4.png")

  i2 <- c(2:4,1) ; nodes2 <- nodes[i2] ; m2 <- m[i2] ; S2 <- S[i2,i2]
  GaussCond(m2,S2,7)

  n <- 5                       ; nodes <- LETTERS[1:n]
  m <- rep(0,n)                ; Vcond <- rep(1,n)
  R <- matrix(0,nrow=n,ncol=n) ; for(j in 2:n){ R[j,j-1] <- 0.5}
  precMatrix( nodes, Vcond, R )

  a  <- 4 ; b  <- -0.22
  mx <- 8 ; Vx <- 1 ; my <- a + b * mx ; m <- c( mx, my ) ; Ve <- 1.e-6
  S  <- matrix( c(  Vx, b  *Vx, b*Vx, b^2*Vx + Ve), nrow=2 )
# Sampling from a multivariate Gaussian that mimics a near-deterministic linear relationship."}
  n <- 100 ; samplexy <- rmvnorm( n, m, S )
  plot( samplexy[,1], samplexy[,2], xlab="x", ylab="y" )
  VR <- VcondR( solve(S) ) ; Vcond <- signif(VR$Vcond,2) ; R <- VR$R

  a  <-  2 ; b1 <- 3 ; b2 <- 5
  mx <- 10 ; Vx <- 1 ; my <-  1 ; Vy <- 1 ; mz <- a + b1 * mx + b2 * my
  m  <- c( mx, my, mz ) ; Ve <- 10
  S  <- matrix( c(    Vx,     0, b1  *Vx,
                       0,    Vy,           b2*  Vy,
                   b1*Vx, b2*Vy, b1^2*Vx + b2^2*Vy + Ve ), nrow=3 )
  n         <- 1e5
  samplexyz <- rmvnorm( n, m, S )
  samplex   <- samplexyz[,1]
  sampley   <- samplexyz[,2]
  samplez   <- samplexyz[,3]
  lm( samplez ~ samplex + sampley )

  n     <- 4 ; nodes <- paste0( "s", 1:4 ) ; m <- rep( 2, n ) ; phi <- 1.9115
  dist  <- matrix( c( 0, 5, 4, 3,
                      5, 0, 3, 4,
                      4, 3, 0, 5,
                      3, 4, 5, 0 ), nrow=4 )
  S     <- 1.47 * exp(-dist/phi)
  VR    <- VcondR( solve(S) )
  Vcond <- signif(VR$Vcond,2) ; R <- signif(VR$R,2)

# DAG for kriging.
  DAGkriging <- grViz( "digraph{ graph[layout = neato]
    node[]
      s1[ label='@@1', pos='0,3!' ] ; s2[ label='@@2', pos='4,0!' ]
      s3[ label='@@3', pos='4,3!' ] ; s4[ label='@@4', pos='0,0!' ]
    edge[]
      s1 -> s2 [ label='@@5-1' ] ; s1 -> s3 [ label='@@5-2' ]
      s1 -> s4 [ label='@@5-3' ] ; s2 -> s3 [ label='@@6-1' ]
      s2 -> s4 [ label='@@6-2' ] ; s3 -> s4 [ label='@@7'   ] }
    [1]: paste0('s1','\\nm=',m[1],',V=',Vcond[1])
    [2]: paste0('s2','\\nm=',m[2],',V=',Vcond[2])
    [3]: paste0('s3','\\nm=',m[3],',V=',Vcond[3])
    [4]: paste0('s4','\\nm=',m[4],',V=',Vcond[4])
    [5]: c( paste0('r=',R[2,1]), paste0('r=',R[3,1]), paste0('r=',R[4,1]) )
    [6]: c( paste0('r=',R[3,2]), paste0('r=',R[4,2]) )
    [7]:    paste0('r=',R[4,3])
  ")
  export_svg(DAGkriging) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGkriging.png")

  i2 <- c(4,1:3) ; nodes2 <- nodes[i2] ; m2 <- m[i2] ; S2 <- S[i2,i2]
  y  <- c(0.7,3.1,2.2) ; GaussCond(m2,S2,y)

# Exercise 2. Covariance matrices (after MacKay).
  M1 = matrix( c(9,  3, 1,  3,  9,  3, 1,  3, 9), nrow=3 )
  M2 = matrix( c(8, -3, 1, -3,  9, -3, 1, -3, 8), nrow=3 )
  M3 = matrix( c(9,  3, 0,  3,  9,  3, 0,  3, 9), nrow=3 )
  M4 = matrix( c(9, -3, 0, -3, 10, -3, 0, -3, 9), nrow=3 )
  VcondR( solve(M1) )
  VcondR( solve(M2) )
  VcondR( solve(M3) )
  VcondR( solve(M4) )
  
# Exercise 3. Simplest hierarchical model.
  Age   <- rbinom( 20, 1, 0.5           )
  Smoke <- rbinom( 20, 1, 0.2+0.6*Age   )
  Run   <- rbinom( 20, 1, 0.4-0.3*Smoke )
  