## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 10. After the Calibration: Interpretation, Reporting, Visualisation

  library(denstrip)
  library(vioplot)

# Four ways of plotting the same distribution. Which one do you prefer?
  par( mfrow=c(1,4), mar=c(3,2,3,2) )
  
  n <- 1e3 ; x <- seq(-4, 4, length=n) ; xsample <- rnorm(n)
  
  dens <- dnorm( x )
  plot(x, xlim=c(-4, 4), ylim=c(-2, 2), xlab="", ylab="", type="n", yaxt="n")
  denstrip(x, dens, at=0, width=1, horiz=T )
  mtext("Density strip")
  
  vioplot( xsample, horizontal=T, drawRect=F, yaxt="n", at=0,
           xaxt="n", ylim=c(-5,5))
  xtick<-seq(-4, 4, by=2)
  axis(side=1, at=xtick )
  mtext( "Violin plot" )
  
  boxplot( xsample, ylim=c(-4,4), horizontal=T, boxwex = 0.5 )
  mtext( "Boxplot" )
  
  hist( xsample, freq=F, main="", xlim=c(-4,4) ) ; box()
  mtext( "Histogram" )
