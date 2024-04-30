## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
### Chapter 12. Discrepancy

  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)

# Error types. Solid line is unknown reality $z$, dashed line is model output
# $f(x,\\theta)$, blue dots are data $y$. Model error
# $\\epsilon_{\\textrm{model}}$ depends on  model structure $f$, conditions $x$
# and model parameterisation $\\theta$. If $\\theta$ is optimal, then the
# remaining model error is called the *discrepancy*; it is due to erroneous
# model structure. Because $z$ is never fully known, we always remain uncertain
# about $\\epsilon_{\\textrm{model}}$ and $\\epsilon_y$.
  set.seed(1)

  n <- 100
  x <- seq( 0, 3, length=n ) ; xobs <- seq( 0.4, 2.8, by = 0.4 )
  z <- function(x) { 10 * (1 - exp(-2*x)) }
  f <- function(x) { 4 * x }
  y <- rnorm( xobs, z(xobs)*0.9, 2 )
  y[1] <- z(xobs[1])+2 ; y[2] <- z(xobs[2])+1
  
  plot   ( x, z(x), type="l", frame.plot=F, ylab="z, f, y", xlim=c(0,3.8), ylim=c(0,13), lwd=2 )
  polygon( c(x[ 1: 83],0), c(z(x[1:83])-0.1,0), col=gray(0.9), border=NA )
  polygon( c(x[83:100],x[100],x[83]), c(z(x[83:100])+0.1,f(x[100]),f(x[83])),
           col=gray(0.9), border=NA )
  points ( x, f(x), type="l", lty=2, lwd=2 )
  points ( xobs, y, col="blue", pch=16 )
  arrows ( 1, f(1), 1, z(1), code=3, length=0.1, lty=1, lwd=2 )
  arrows ( xobs[1], y[1], xobs[1], z(xobs[1]), code=3, length=0.1, lwd=2, col="blue" )
  arrows ( xobs[4], f(xobs[4]), xobs[4], y[4], code=3, length=0.1, lwd=2, col="red" )
  text(3.15 ,10  , "z"                                , col="black", cex=1.3)
  text(3.3 ,12.5, expression(paste("f(x,",theta,")")), col="black", cex=1.3)
  text(0.25, 6.5, expression(epsilon[y]       )      , col="blue" , cex=1.3)
  text(0.75, 5  , expression(epsilon[model]   )      , col="black", cex=1.3)
  text(1.25,10.5, expression(epsilon[mismatch])      , col="red"  , cex=1.3)

# Conceptual diagram for model discrepancy, showing model-specific biases and an
# overall bias shared by the whole cluster. Simplified after Chandler (2013).
  Chandler <- grViz("digraph neato {
    graph[layout=neato]
    node [shape=circle, style=filled, color=grey, fillcolor=lightblue]
      z[label='Reality']
    node [shape=point]
      b
    node [shape=circle, fillcolor = orange]
    edge [color = grey, len=4, label= '@@6']
      z -> {b}
    edge [color = grey, len=1]
      b -> M1 [label='@@1']
      b -> M2 [label='@@2']
      b -> M3 [label='@@3']
    edge [color = grey, len=1.5]
      b -> M4 [label='@@4']
      b -> M5 [label='@@5'] }
    [1]: paste0('&epsilon;&#x2081;')
    [2]: paste0('&epsilon;&#x2082;')
    [3]: paste0('&epsilon;&#x2083;')
    [4]: paste0('&epsilon;&#x2084;')
    [5]: paste0('&epsilon;&#x2085;')
    [6]: paste0('&epsilon;&#x2080;')
  ")
  export_svg(Chandler) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("Chandler.png")
