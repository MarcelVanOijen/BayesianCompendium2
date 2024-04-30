## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 5. Deriving the Posterior Distribution {#ChPosterior}

  m0 <- 0               ; V0 <- 1
  y  <- c(1.0,0.9,1.1)  ; ny <- length(y)
  my <- mean(y)         ; Vy <- 1
  k  <- (1/V0) / (1/V0 + ny*Vy)
  m1 <- k*m0 + (1-k)*my ; V1 <- 1 / (1/V0 + ny/Vy)

  Lmax       <- prod( dnorm( y, my, sqrt(Vy) ) )
  n          <- 1e5 ; samplePrior <- rnorm( n, m0, sqrt(V0) )
  L          <- function(b){ prod( dnorm( y, b, sqrt(Vy) ) ) }
  Lrelative  <- sapply(samplePrior,L) / Lmax
  iAccept    <- which( runif(n) < Lrelative )
  samplePost <- samplePrior[ iAccept ]
  m1.n       <- mean( samplePost ) ; V1.n <- var( samplePost )

# Sampling from the posterior distribution of this chapter by means of Accept-Reject."}
  hist( samplePost, xlab="", ylab="", main=paste0( "A-R sample from posterior\n(Mean = ",
                    signif(m1.n,3) , "; Var = ", signif(V1.n,3), ")" ) )
