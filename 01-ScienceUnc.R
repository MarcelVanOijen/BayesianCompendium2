## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 1. Science and Uncertainty

  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)

# The four model types. The deterministic models can be made stochastic by
# embedding them in a statistical model.
  DAGmodels <- grViz( "digraph{ graph[]
    node[ shape=box ]
      A [label='@@1'] ; B1[label='@@2'] ; B2[label='@@3']
      C1[label='@@4'] ; C2[label='@@5'] ; C3[label='@@4'] ; C4[label='@@5']
    edge[]
      A -> B1 ; A -> B2 ; B1 -> C1 ; B1 -> C2 ; B2 -> C3 ; B2 -> C4 }
    [1]: 'Models'
    [2]: 'Process-based'
    [3]: 'Empirical'
    [4]: 'Deterministic'
    [5]: 'Stochastic'
  ")
  export_svg(DAGmodels) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGmodels.png")

# Exercise 1. Uncertainty propagation.
  fa <- function(x,b) {x + b} ; fb <- function(x,b) {x - b}
  fc <- function(x,b) {x * b} ; fd <- function(x,b) {exp(-(x + b)**2)}
  n  <- 1e4 ; x <- rnorm(n,1,1) ; b <- rnorm(n,2,1)
  par( mfrow=c(1,4) ) # This command prepares 1 row of 4 plots
  hist( fa(x,b) ) ; hist( fb(x,b) ) ; hist( fc(x,b) ) ; hist(fd(x,b) )
  