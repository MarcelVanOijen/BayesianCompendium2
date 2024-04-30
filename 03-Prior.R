## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Ch. 3. Assigning a Prior Distribution

  library(FD)
  Dp <- 1:6 ; meA <- FD::maxent(3.5, Dp)
  pA <- meA$prob ; SA <- meA$entropy
  SA <- -sum(pA * log(pA)) # = meA$entropy

  meB      <- maxent(4.5, Dp)                     ; pB <- meB$prob
  meC      <- maxent( 13, Dp^2 )                  ; pC <- meC$prob
  meD      <- maxent( c(3.5,16), rbind(Dp,Dp^2) ) ; pD <- meD$prob
  meE      <- maxent( c(3.5,13), rbind(Dp,Dp^2) ) ; pE <- meE$prob
  entropyB <- meB$entropy ; SB <- -sum(pB * log(pB))
  entropyC <- meC$entropy ; SC <- -sum(pC * log(pC))
  entropyD <- meD$entropy ; SD <- -sum(pD * log(pD))
  entropyE <- meE$entropy ; SE <- -sum(pE * log(pE))

# Distributions selected using MaxEnt. Entropies between brackets.
# All distributions are on the same domain (1,..,6) but with different
# constraints for mean and mean square: see text.
  par(mfrow=c(2,3),mar=c(2,2,3,1))
  barplot( pA, main=paste0("A\n(S=",signif(SA,4),")"), names=1:6 )
  barplot( pB, main=paste0("B\n(S=",signif(SB,4),")"), names=1:6 )
  barplot( pC, main=paste0("C\n(S=",signif(SC,4),")"), names=1:6 )
  barplot( pD, main=paste0("D\n(S=",signif(SD,4),")"), names=1:6 )
  barplot( pE, main=paste0("E\n(S=",signif(SE,4),")"), names=1:6 )

  DpF      <- seq(1,6,length.out=51)
  meF      <- maxent( c(3.5,13), rbind(DpF,DpF^2) ) ; pF <- meF$prob
  entropyF <- meF$entropy ; SF  <- -sum(pF * log(pF))
  DpG      <- 1:10
  meG      <- maxent( c(6,100.5), rbind( DpG, 100+(DpG-5.5)^3 ) ) ; pG <- meG$prob
  entropyG <- meG$entropy
  SG       <- -sum(pG * log(pG))

# More MaxEnt distributions. Domains cover 51 and 10 possible values.
# Entropies between brackets.
  par(mfrow=c(1,2),mar=c(3,2,3,1))
  barplot( pF, main=paste0("F\n(S=",signif(SF,4),")"),
           names=1 + (0:50)/10, cex.names=0.5 )
  barplot( pG, main=paste0("G\n(S=",signif(SG,4),")"),
           names=1:10, cex.names=0.5 )
  
# Exercise 1. Family matters (prior).
  p.b <- 1 / sum( 1/(1:6) )
  p.c <- sum( 1/(1:6)^2 ) / sum( 1/(1:6) )
