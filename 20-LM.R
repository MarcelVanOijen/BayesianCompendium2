## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 20. Linear Modelling: LM, GLM, GAM and Mixed Models

  library(arm)
  library(brms)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(mgcv)
  library(rstanarm)
  library(rsvg)

# Generalisation of the basic linear model (LM) to generalised linear and
# additive models (GLM, GAM) and mixed models (LMM, GLMM, GAMM).
  DAGlinearmodels <- grViz( "digraph{ graph[layout = neato]
    node[ shape=box, width=2 ]
      s1[ label='@@1', pos='0,3!', fontsize=40 ]
      s2[ label='@@2', pos='4,3!', fontsize=40 ]
      s3[ label='@@3', pos='8,3!', fontsize=40 ]
      s4[ label='@@4', pos='0,0!', fontsize=40 ]
      s5[ label='@@5', pos='4,0!', fontsize=40 ]
      s6[ label='@@6', pos='8,0!', fontsize=40 ]
    edge[]
      s1 -> s2 -> s3 ; s1 -> s4 ; s2 -> s5 ; s3 -> s6
      s4 -> s5 -> s6 }
    [1]: paste0('LM'  )
    [2]: paste0('GLM' )
    [3]: paste0('GAM' )
    [4]: paste0('LMM' )
    [5]: paste0('GLMM')
    [6]: paste0('GAMM')
  ")
  export_svg(DAGlinearmodels) %>% charToRaw() %>% rsvg() %>%
    png::writePNG("DAGlinearmodels.png")

  set.seed(1)
  x <- 0:9 ; y <- exp(x/4) + sin(x) + rnorm(10)
  xy <- data.frame(x,y)
  lm.1 <- lm      ( y~x )
  lm.2 <- glm     ( y~x, family=gaussian(link=identity) )
  lm.3 <- bayesglm( y~x, family=gaussian(link=identity) )
  lm.4 <- brm     ( y~x, family=gaussian(link=identity), data=xy )
  lm.5 <- stan_glm( y~x, family=gaussian(link=identity), data=xy )
  lm.6 <- gam     ( y~x, family=gaussian(link=identity) )

  set.seed(1)
  i0    <- c(1,3,4:6) ; i1 <- setdiff(1:10,i0)
  group <- rep(0,length(x)) ; group[i1] <- 1
  xgy   <- data.frame(x,group,y)
  lmm.1 <- brm     ( y~x+(1|group), data=xgy )
  glm.1 <- stan_glm( y~x, family=gaussian(link=log), data=xy )
  gam.1 <- gam     ( y~s(x) )

# LM, LMM with varying intercept and common slope, GLM with log link-function,
# GAM with thin plate splines. Each model is applied to the same dataset
# ($x$,$y$). In brackets: software used.
  par(mfrow=c(1,4), mar=c(3,2,3,1) )
  plot( x, y, col="blue", ylim=c(0,10), main="LM\n(arm)" )
  pred  <- predict(object=lm.3, type="response")
  lines( x, pred )
  plot( x, y, col="blue", ylim=c(0,10), main="LMM\n(brm)", yaxt="n" )
  pred   <- predict(object=lmm.1, type="response")
  points( x[i1], y[i1], pch=16 )
  lines ( x[i0], pred[i0,1] ) ; lines ( x[i1], pred[i1,1] )
  legend( "topleft", legend=c("group=1","group=0"), pch=c(16,1), col="darkblue" )
  plot( x, y, col="blue", ylim=c(0,10), main="GLM\n(rstanarm)", yaxt="n" )
  pred  <- predict(object=glm.1, type="response")
  lines( x, pred )
  plot( x, y, col="blue", ylim=c(0,10), main="GAM\n(mgcv)", yaxt="n" )
  pred  <- predict(object=gam.1, type="response")
  lines( x, pred )
