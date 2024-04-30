## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 19. Bayesian Decision Theory

  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)

# A graphical model for Bayesian decision theory.
  DAGBDT <- grViz( "digraph {
    graph [rankdir=BT]
  
    node [shape=circle, fixedsize = true, width = 1.2]
    a [label='@@1'] ; e [label='@@2'] ; C [label='@@3']
    t [label='@@4'] ; u [label='@@5'] ; x [label='@@6']
    z [label='@@7'] ; B [label='@@8']
    edge [minlen=1]
    nodesep = 0.7
    {rank = same; t, x, a}
    {rank = same; z, e}
    {rank = same; C, B}
    a      -> C -> u
    t -> z -> B -> u
    a -> z
    x -> z
    e -> z
    C [fillcolor=lightgray,style=filled]
    a [shape=box,fillcolor=lightgray,style=filled,width=1.6]
    u [shape=ellipse,fillcolor=lightgray,style=filled,width=2] }
    [1]: paste0( 'Action'         , '\\na'         )
    [2]: paste0( '\\nError'       , '\\n&epsilon;' )
    [3]: paste0( 'Cost'                            )
    [4]: paste0( '\\nParameters'  , '\\n&theta;'   )
    [5]: paste0( 'Utility'        , '\\nu'         )
    [6]: paste0( '\\nEnvironment' , '\\nx'         )
    [7]: paste0( '\\nSystem'      , '\\nz'         )
    [8]: paste0( 'Benefit'                         )
  ")
  export_svg(DAGBDT) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGBDT.png")

# Example of linear BDT represented by a DAG. Nodes show the values of the prior
# mean and conditional variance. Edges are labelled with the regression
# coefficient.
  nodes <-         c( "IRRIG", "RAIN", "YIELD", "COST", "BENEFIT", "UTILITY")
  m     <-         c(  0     ,  500  ,   5    ,    0  ,  150     ,  150     )
  Vcond <-         c(  0.01  ,  1e4  ,   1    ,  100  ,  100     ,  100     )
  R     <- matrix( c(  0     ,  0    ,   0    ,    0  ,    0     ,    0 ,
                       0     ,  0    ,   0    ,    0  ,    0     ,    0 ,
                       0.01  ,  0.01 ,   0    ,    0  ,    0     ,    0 ,
                       0.5   ,  0    ,   0    ,    0  ,    0     ,    0 ,
                       0     ,  0    ,  30    ,    0  ,    0     ,    0 ,
                       0     ,  0    ,   0    ,   -1  ,    1     ,    0     ),
                   byrow=T, nrow=6 )
  DAGBDTlinear <- grViz( "digraph{
    graph[rankdir=BT]
    node[]
      IRRIG[label='@@1'] ; RAIN   [label='@@2'] ; YIELD  [label='@@3']
      COST [label='@@4'] ; BENEFIT[label='@@5'] ; UTILITY[label='@@6']
    {rank = same; IRRIG, RAIN}
    {rank = same; YIELD}
    {rank = same; COST, BENEFIT}
    {rank = same; UTILITY}
    edge[]
      IRRIG -> YIELD  [label='@@7-1'] ; RAIN    -> YIELD  [label='@@7-2']
      YIELD -> BENEFIT[label='@@8-1'] ; IRRIG   -> COST   [label='@@8-2']
      COST  -> UTILITY[label='@@9-1'] ; BENEFIT -> UTILITY[label='@@9-2']
    IRRIG  [shape=box,fillcolor=lightgray,style=filled]
    COST   [          fillcolor=lightgray,style=filled]
    UTILITY[          fillcolor=lightgray,style=filled] }
    [1]: paste0( 'IRRIG'  , '\\nm=', m[1], ',V=', Vcond[1] )
    [2]: paste0( 'RAIN'   , '\\nm=', m[2], ',V=', Vcond[2] )
    [3]: paste0( 'YIELD'  , '\\nm=', m[3], ',V=', Vcond[3] )
    [4]: paste0( 'COST'   , '\\nm=', m[4], ',V=', Vcond[4] )
    [5]: paste0( 'BENEFIT', '\\nm=', m[5], ',V=', Vcond[5] )
    [6]: paste0( 'UTILITY', '\\nm=', m[6], ',V=', Vcond[6] )
    [7]: c( paste0('r=',R[3,1]), paste0('r=',R[3,2]) )
    [8]: c( paste0('r=',R[5,3]), paste0('r=',R[4,1]) )
    [9]: c( paste0('r=',R[6,4]), paste0('r=',R[6,5]) )
  ")
  export_svg(DAGBDTlinear) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGBDTlinear.png")

  u <- function( a, x=1, t=1, e=0, ka=0.2, kz=1 ) {
    z       <- t * (1-exp(-a-x)) + e
    cost    <- ka * a
    benefit <- kz * z
    return( benefit - cost ) }
  na    <- 21
  a.seq <- seq( 0, 2, length.out=na )
  u.seq <- u(a.seq)

  set.seed(1) ; np <- 1e3
  x.smp1  <- rnorm( np, 1  , 1   ) ; t.smp1  <- rnorm( np, 1  , 0.5 )
  e.smp1  <- rnorm( np, 0  , 1   )
  ka.smp1 <- runif( np, 0.1, 0.3 ) ; kz.smp1 <- runif( np, 0.5, 1.5 )
  kunc    <- 1.5
  x.smp2  <- 1   + (x.smp1  - 1   ) * kunc
  t.smp2  <- 1   + (t.smp1  - 1   ) * kunc
  e.smp2  <-        e.smp1          * kunc
  ka.smp2 <- 0.2 + (ka.smp1 - 0.2 ) * kunc
  kz.smp2 <- 1   + (kz.smp1 - 1   ) * kunc
  
  u.tbl1  <- u.tbl2 <- NULL
  for(i in 1:np) {
    u.tbl1 <- rbind( u.tbl1,
      u( a.seq, x.smp1[i], t.smp1[i], e.smp1[i], ka.smp1[i], kz.smp1[i] ) )
    u.tbl2 <- rbind( u.tbl2,
      u( a.seq, x.smp2[i], t.smp2[i], e.smp2[i], ka.smp2[i], kz.smp2[i] ) )
  }

  Qlo       <- 0.25 ; Qhi <- 0.75
  uQlo.seq1 <- sapply( 1:na, function(i) { quantile( u.tbl1[,i], probs=Qlo ) } )
  uQlo.seq2 <- sapply( 1:na, function(i) { quantile( u.tbl2[,i], probs=Qlo ) } )
  umn.seq1  <- colMeans(u.tbl1)
  umn.seq2  <- colMeans(u.tbl2)
  uQhi.seq1 <- sapply( 1:na, function(i) { quantile( u.tbl1[,i], probs=Qhi ) } )
  uQhi.seq2 <- sapply( 1:na, function(i) { quantile( u.tbl2[,i], probs=Qhi ) } )
  imaxmnu.1 <- which( umn.seq1 == max(umn.seq1) )
  amaxmnu.1 <- a.seq[imaxmnu.1]
  imaxmnu.2 <- which( umn.seq2 == max(umn.seq2) )
  amaxmnu.2 <- a.seq[imaxmnu.2]

# Utility as a function of environment $x$ for a negative exponential
# performance function and linear cost and benefit functions.
# Left: action $a = 0$; right $a = 2$.
# Both panels are for the same sample (n = 1000) of $x$ and parameters.
  par( mfrow=c(1,2), mar=c(4,2,3,0) )
  u.range <- range( u.tbl1[,1], u.tbl1[,21] )
  plot( x.smp1, u.tbl1[, 1], ylim=u.range, xlab="x", ylab="",
        main ="Utility\n(a=0)" )
  plot( x.smp1, u.tbl1[,na], ylim=u.range, xlab="x", ylab="", yaxt="n",
        main ="Utility\n(a=2)" )

# Expected utility as a function of action $a$ for a negative exponential
# performance function and linear cost and benefit functions.
# Left panel  : no parameter uncertainty.
# Middle panel: default uncertainties (see text).
# Right panel : uncertainties multiplied by 1.5.
# Solid line  : expectation.
# Dashed lines: Q25 and Q75.
  par( mfrow=c(1,3), mar=c(4,2,3,0) )
  u.range <- range( uQlo.seq1, umn.seq1, uQhi.seq1,
                    uQlo.seq2, umn.seq2, uQhi.seq2 )
  imaxu <- which( u.seq == max(u.seq) ) ; amaxu <- a.seq[imaxu]
  plot  ( a.seq, u.seq, xlab="a", ylab="", main="Utility\n(no uncertainty)",
          type="l", ylim=u.range )
  abline( v=amaxu, col="red", lty=2 )
  text( amaxu, 0, labels=paste("a* =",amaxu ), pos=4 )
  plot  ( a.seq, umn.seq1 , xlab="a", ylab="",
          main="Utility\n(low uncertainty)", type="l", ylim=u.range, yaxt="n" )
  points( a.seq, uQlo.seq1, type="l", lty=2 )
  points( a.seq, uQhi.seq1, type="l", lty=2 )
  abline( v=amaxmnu.1, col="red", lty=2 )
  text( amaxmnu.1, 0, labels=paste("a* =",amaxmnu.1), pos=4 )
  plot  ( a.seq, umn.seq2 , xlab="a", ylab="",
          main="Utility\n(high uncertainty)", type="l", ylim=u.range, yaxt="n" )
  points( a.seq, uQlo.seq2, type="l", lty=2 )
  points( a.seq, uQhi.seq2, type="l", lty=2 )
  abline( v=amaxmnu.2, col="red", lty=2 )
  text( amaxmnu.2, 0, labels=paste("a* =",amaxmnu.2), pos=2 )

  PRA <- function( x, z, thr=0 ) {
    n  <- length(z) ; H    <- which(x < thr) ; nH      <- length(H)
    Ez <- mean( z ) ; Ez_H <- mean( z[H] )   ; Ez_notH <- mean( z[-H] )
    pH <- nH / n    ; V    <- Ez_notH - Ez_H ; R       <- Ez_notH - Ez
    return( c(pH=pH,V=V,R=R) ) }
# PRAs on data generated for BDT.
  thr      <- 0
  PRA.BDT1 <- sapply( 1:na, function(i){ PRA( x.smp1+a.seq[i], u.tbl1[,i], thr ) } )
  PRA.BDT2 <- sapply( 1:na, function(i){ PRA( x.smp2+a.seq[i], u.tbl2[,i], thr ) } )
  par( mfrow=c(1,3), mar=c(4,2,3,0) )
  pHlim <- range( c( 0, PRA.BDT1["pH",], PRA.BDT2["pH",] )          )
  Vlim  <- range( c( 0, PRA.BDT1["V" ,], PRA.BDT2["V" ,] ), na.rm=T )
  plot  ( a.seq, PRA.BDT1["pH",], type="l", ylim=pHlim,
          main="p[H]", xlab="a", ylab="")
  points( a.seq, PRA.BDT2["pH",], type="l", lty=2 )
  legend( "topright", legend=c("High unc", "Low unc"), lty=2:1, cex=0.75 )
  plot  ( a.seq, PRA.BDT1["V",], col="black", type="l", ylim=Vlim,
          main="V", xlab="a", ylab="" )
  points( a.seq, PRA.BDT2["V",], col="black", type="l", lty=2 )
  legend( "topright", legend=c("High unc", "Low unc"), lty=2:1, cex=0.75 )
  plot  ( a.seq, PRA.BDT1["R",], col="black", type="l", ylim=Vlim,
          main="R", xlab="a", ylab="" )
  points( a.seq, PRA.BDT2["R",], col="black", type="l", lty=2 )
  legend( "topright", legend=c("High unc", "Low unc"), lty=2:1, cex=0.75 )

# Impact of action $a$ on corrected risk $R_c$ calculated from data generated
# for BDT. $R_c$ reaches an uncertainty-dependent minimum at $a*$.
  fz <- function(a,x,t) { t * (1-exp(-a-x)) }
  z.tbl1 <- z.tbl2 <- NULL
  for(i in 1:np) {
    z.tbl1 <- rbind( z.tbl1, fz(a.seq,x.smp1[i],t.smp1[i]) + e.smp1[i] )
    z.tbl2 <- rbind( z.tbl2, fz(a.seq,x.smp2[i],t.smp2[i]) + e.smp2[i] )
  }
  h.lst1 <- h.lst2 <- vector( "list", na )
  for(i in 1:na) {
    h.lst1[[i]] <- which( (a.seq[i]+x.smp1) > thr )
    h.lst2[[i]] <- which( (a.seq[i]+x.smp2) > thr )
  }
  benefits.mnabove.seq1 <- benefits.mnabove.seq2 <- rep(NA,na)
  costs.mnabove.seq1    <- costs.mnabove.seq2    <- rep(NA,na)
  for(i in 1:na) {
    benefits.mnabove.seq1[i] <- mean( (kz.smp1*z.tbl1[,i])[ h.lst1[[i]] ] )
    benefits.mnabove.seq2[i] <- mean( (kz.smp2*z.tbl2[,i])[ h.lst2[[i]] ] )
    costs.mnabove.seq1   [i] <- a.seq[i] * mean( ka.smp1  [ h.lst1[[i]] ] )
    costs.mnabove.seq2   [i] <- a.seq[i] * mean( ka.smp2  [ h.lst2[[i]] ] )
  }
  Rc.seq1  <- PRA.BDT1["R",] + costs.mnabove.seq1 - benefits.mnabove.seq1
  Rc.seq2  <- PRA.BDT2["R",] + costs.mnabove.seq2 - benefits.mnabove.seq2
  imaxRc.1 <- which( Rc.seq1 == min(Rc.seq1,na.rm=T) ) ; amaxRc.1 <- a.seq[imaxRc.1]
  imaxRc.2 <- which( Rc.seq2 == min(Rc.seq2,na.rm=T) ) ; amaxRc.2 <- a.seq[imaxRc.2]
  par( mfrow=c(1,1) )
  Rclim  <- range( 0, Rc.seq1, Rc.seq2, na.rm=T )
  plot  ( a.seq, Rc.seq1, type="l",
          main="Rc", xlab="a", ylab="", ylim=Rclim )
  points( a.seq, Rc.seq2, type="l", lty=2 )
  abline( v=amaxRc.1              , lty=1 )
  abline( v=amaxRc.2              , lty=2 )
  text( amaxRc.1, 0.05+min(Rc.seq1,na.rm=T), labels=paste("a* =",amaxRc.1 ), pos=4 )
  text( amaxRc.2, 0.05+min(Rc.seq2,na.rm=T), labels=paste("a* =",amaxRc.2 ), pos=4 )
  legend( "topright", legend=c("High unc", "Low unc"), lty=2:1, cex=0.75 )
