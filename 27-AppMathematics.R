## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Appendix B: Mathematics for Modellers

  p <- function(b,m=1,s=1) {exp( -0.5 * ((b-m)/s)^2 ) / (s*sqrt(2*pi))}
  b <- seq(-10,10,by=0.1) ; plot( b, p(b), type="l" )
