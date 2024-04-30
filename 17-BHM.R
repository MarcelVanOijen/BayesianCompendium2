## M. van Oijen (2024). Bayesian Compendium, 2nd edition.
## Chapter 17. Bayesian Hierarchical Modelling {#ChBHM}

  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(ggplot2)
  library(gridExtra)
  library(rsvg)
  library(rjags)

# Three statistical models for n observations on J different grassland
# cultivars. From left to right: Model A is a non-hierarchical model with
# universal slope ('a') and intercept ('b'); Model B is also non-hierarchical
# but has cultivar-specific parameters a and b; Model C is a Bayesian
# hierarchical model (BHM) with cultivar-specific parameters and global
# hyperparameters. We use plate-notation where the 'n' and 'J' in bottom-right
# corners indicate the multiplicity of data and cultivars.
# Model A.
  DAGM1 <- grViz("digraph{ graph[rankdir=BT, ranksep=1, nodesep=0.5]
      node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        a b
      subgraph cluster_1 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='@@1']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        y
      }
      edge[arrowhead=vee, arrowsize=1.25]
        a -> y ; b -> y
      }
      [1]: '.               n'
    ")
  export_svg(DAGM1) %>% charToRaw() %>% rsvg() %>% png::writePNG("DAGM1.png")
# Model B.
  DAGM2 <- grViz("digraph{ graph[rankdir=BT, ranksep=1, nodesep=0.5]
      node[fillcolor=White, margin=0.1, shape=circle, style=filled]
      subgraph cluster_1 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='@@1']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        b a
      }
      subgraph cluster_2 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='@@2']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        y
      }
      edge[arrowhead=vee, arrowsize=1.25]
        a -> y ; b -> y
      }
      [1]: '.                              J'
      [2]: '.               n'
    ")
  export_svg(DAGM2) %>% charToRaw() %>% rsvg() %>% png::writePNG("DAGM2.png")
# Model C.
  DAGM3 <- grViz("digraph{ graph[rankdir=BT, ranksep=1, nodesep=0.5]
      node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        # xia[label = <&xi;<FONT POINT-SIZE='10'><SUB>a</SUB></FONT>>]
        xia[label = <&xi;<SUB>a</SUB>>]
        xib[label = <&xi;<SUB>b</SUB>>]
      subgraph cluster_1 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='@@1']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        a b
      }
      subgraph cluster_2 {
        graph[shape=rectangle, style=rounded, bgcolor=LemonChiffon, label='@@2']
        node[fillcolor=White, margin=0.1, shape=circle, style=filled]
        y
      }
      edge[arrowhead=vee, arrowsize=1.25]
        xia -> a ; xib -> b ; a -> y ; b -> y
      }
      [1]: '.                              J'
      [2]: '.               n'
    ")
  export_svg(DAGM3) %>% charToRaw() %>% rsvg() %>% png::writePNG("DAGM3.png")

  t  <- c(141, 148, 154, 162, 169, 132, 139, 146, 153, 160, 167, 173, 131, 138, 145, 152, 159, 166, 173, 131, 138, 145, 152, 159, 166, 173, 131, 138, 145, 152, 159, 166, 173, 131, 138, 145, 152, 159, 166, 173, 139, 146, 152, 159, 167, 173, 180, 187, 195, 201, 207, 216, 221, 229, 236, 243, 257, 273, 125, 130, 136, 144, 150, 156, 164, 171, 178, 185, 192, 199, 206, 214, 220, 227, 234, 241, 164, 167, 195, 209, 223, 162, 169, 197, 211, 225, 153, 161, 167, 209, 153, 161, 167, 209, 153, 161, 167, 209, 158, 165, 172, 214, 158, 165, 172, 214, 158, 165, 172, 214, 157, 164, 171, 219, 157, 164, 171, 219, 157, 164, 171, 219, 168, 172, 174, 179, 208, 221, 236, 168, 172, 174, 179, 208, 221, 236, 165, 167, 178, 181, 207, 223, 237, 165, 167, 178, 181, 207, 223, 237, 164, 170, 172, 177, 206, 221, 235, 164, 170, 172, 177, 206, 221, 235, 124, 131, 139, 150, 151, 165, 173, 185, 195, 206, 216, 124, 131, 139, 150, 158, 165, 173, 174, 185, 195, 206, 216, 227, 238, 115, 130, 141, 150, 164, 165, 176, 186, 197, 211, 225, 239, 115, 130, 141, 150, 164, 176, 186, 192, 197, 211, 225, 239, 109, 119, 128, 135, 143, 150, 151, 162, 169, 178, 186, 196, 203, 210, 224, 238, 109, 119, 128, 135, 143, 150, 162, 169, 178, 186, 196, 203, 210, 224, 238, 160, 166, 171, 205, 214, 221, 235, 166, 176, 179, 206, 218, 235, 246)
  y  <- c(235.7, 464.1, 681.3, 748.8, 694.6, 84.2, 189.6, 225.4, 365.8, 559, 584, 905, 75, 69.7, 105.2, 126.9, 174.5, 197, 275, 84.7, 117.3, 225.3, 310.5, 454.1, 567, 739, 102.4, 176.1, 257, 367, 503.9, 630, 832, 96.4, 149.2, 264.6, 385.1, 483.8, 640, 865, 178.2, 213.7, 244.1, 261.5, 410.6, 419.6, 566.6, 657.7, 663.8, 515.8, 136.9, 127.8, 119.3, 155.9, 131.7, 188.9, 172.5, 177.3, 173.6, 95.4, 152.6, 224.2, 129.3, 275.7, 321.4, 282, 328.6, 268.1, 51.3, 124.4, 166.1, 148.1, 68.3, 107.2, 87.7, 222.5, 466.3, 559.9, 317.6, 422.7, 483.4, 510.6, 655.6, 237.8, 415.9, 509.8, 324.6, 472.6, 770.6, 226.6, 448.6, 666.6, 787.6, 191.6, 324, 515, 911, 192, 314.6, 384.6, 417.6, 162.6, 471.6, 685.6, 725.6, 205.6, 542.6, 697.6, 903.6, 330.6, 243.6, 393.6, 454.6, 162.6, 388.6, 662.6, 674.6, 130.6, 475.6, 660.6, 834.6, 187.6, 174.8, 226.8, 288.5, 378.6, 272, 355.2, 433.9, 212.4, 273.9, 298.1, 418.1, 280.8, 382.4, 473.7, 272.6, 294.8, 485.2, 590.2, 204.4, 397.4, 498.9, 280.1, 308.4, 498.4, 534.2, 222.2, 415.2, 500.8, 181.7, 300.7, 408, 533.7, 218.1, 292.5, 373.1, 173.9, 297.4, 420, 526, 182.5, 350.6, 443.4, 132.3, 278.9, 434.1, 686.4, 90.9, 70.7, 95.8, 215.6, 362.1, 588.8, 768.1, 132.3, 278.9, 434.1, 686.4, 844, 1151.8, 1172, 202, 82.9, 147.3, 352.2, 518.8, 624.3, 767.4, 36.2, 135.1, 257.3, 553.3, 702.2, 71.9, 29.7, 107.9, 309, 546.5, 681.5, 626.7, 36.2, 135.1, 257.3, 553.3, 702.2, 1237.4, 1270.1, 95, 75.8, 233.3, 354.6, 373.6, 20, 83.4, 153.6, 280.3, 470.7, 567.7, 100.6, 15.4, 55.5, 143.7, 234.1, 308.2, 416.8, 526.1, 726.1, 580.3, 20, 83.4, 153.6, 280.3, 470.7, 567.7, 762.1, 757.3, 813.9, 39, 118.6, 190.3, 255.4, 470.2, 559.5, 329.9, 437.5, 541.5, 221.5, 376.3, 409.9, 594.1, 251.4, 401.5, 439.9, 208.3, 378.6, 553.9, 617.3)
  cv <- c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)
  s  <- c(28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27)
  ny <- length(y) ; ncv <- length(unique(cv))
# summ1 <- aggregate( y, list(cv), summary )
# summtable <- cbind(summ1[,1],summ1[,2]) ; colnames(summtable)[1] <- "cv"

# Grassland growth data used for Bayesian calibration of three models.
# 't' = time (days), 'cv' = cultivar.
  df_DM <- data.frame(t,y,cv) ; df_DM$cv <- as.factor(df_DM$cv)
  p1        <- ggplot(df_DM, aes(t, y))
  settings1 <- list( geom_point(aes(colour = cv)), xlim(100, 300), ylab(""),
                     theme(legend.position="none"), ggtitle("Biomass (g m-2)") )
  p1a <- p1  + settings1
  p1b <- p1a + facet_wrap(~ cv) + theme(legend.position="right") +
               theme(axis.text.x=element_blank())
  grid.arrange(p1a,p1b,ncol=2)

# ModelA
  ModelA <- " model {
    for (i in 1:ny) {
      y[i] ~ dnorm( a + b * t[i], tau.y[i] )
      tau.y[i] <- pow( sigma.y[i], -2 ) }
    a ~ dnorm( 5e2, 1e-6 ) ; b ~ dnorm( 1, 1e-2 ) } "
  writeLines( ModelA, con="ModelA.txt" )
  data.A <- list ( ny=ny, y=y, t=t, sigma.y=rep(1e2,ny) )
  ModelA <- jags.model( "ModelA.txt", data=data.A, n.chains=3, n.adapt=1e4 )
  update( ModelA, n.iter=1e4 )
  ModelA.coda <- coda.samples( ModelA, var=c("a","b"), n.iter=1e4 )
# Post-processing output from JAGS (ModelA)
# summary( ModelA.coda )
  ModelA.mcmc <- as.matrix( ModelA.coda )
  i.a <- which(colnames(ModelA.mcmc)=="a") ; i.b <- which(colnames(ModelA.mcmc)=="b")
  ModelA.a <- ModelA.mcmc[ , i.a ] ; ModelA.b <- ModelA.mcmc[ , i.b ]

# ModelB
  ModelB <- " model {
    for (i in 1:ny) {
      y[i]      ~ dnorm( a[cv[i]] + b[cv[i]] * t[i], tau.y[i] )
      tau.y[i] <- pow( sigma.y[i], -2 ) }
    for (c in 1:ncv) { a[c] ~ dnorm( 5e2, 1e-6 ) ; b[c] ~ dnorm( 1, 1e-2 ) } } "
  writeLines( ModelB, con="ModelB.txt" )
  data.B <- list ( ny=ny, ncv=ncv, y=y, t=t, cv=cv, sigma.y=rep(1e2,ny) )
  ModelB <- jags.model( "ModelB.txt", data=data.B, n.chains=3, n.adapt=1e4 )
  update( ModelB, n.iter=1e4 )
  ModelB.coda <- coda.samples( ModelB, var=c("a","b"), n.iter=1e4 )
# Post-processing output from JAGS (ModelB)
# summary( ModelB.coda )
  ModelB.mcmc <- as.matrix( ModelB.coda )
  i.a      <- sapply( 1:ncv, function(i) {
    which( colnames(ModelB.mcmc)==paste("a[",i,"]",sep="") ) } )
  i.b      <- sapply( 1:ncv, function(i) {
    which( colnames(ModelB.mcmc)==paste("b[",i,"]",sep="") ) } )
  ModelB.a <- ModelB.mcmc[ , i.a ] ; ModelB.b <- ModelB.mcmc[ , i.b ]

# ModelC
  ModelC <- " model {
    for (i in 1:ny) {
      y[i]      ~ dnorm( a[cv[i]] + b[cv[i]] * t[i], tau.y[i] )
      tau.y[i] <- pow( sigma.y[i], -2 ) }
    for (c in 1:ncv) { a[c] ~ dnorm( mu.a, tau.a ) ; b[c] ~ dnorm( mu.b, tau.b ) }
    mu.a ~ dnorm(5e2,1e-6) ; tau.a <- pow(sigma.a,-2) ; sigma.a ~ dunif(0,1e3)
    mu.b ~ dnorm(1,1e-6)   ; tau.b <- pow(sigma.b,-2) ; sigma.b ~ dunif(0,2) } "
  writeLines( ModelC, con="ModelC.txt" )
  data.C <- list ( ny=ny, ncv=ncv, y=y, t=t, cv=cv, sigma.y=rep(1e2,ny) )
  ModelC.params <- c("a","b","mu.a","sigma.a","mu.b","sigma.b")
  ModelC <- jags.model( "ModelC.txt", data=data.C, n.chains=3, n.adapt=1e4 )
  update( ModelC, n.iter=1e4 )
  ModelC.coda <- coda.samples( ModelC, var=ModelC.params, n.iter=1e4 )
# Post-processing output from JAGS (ModelC)
# summary( ModelC.coda )
  ModelC.mcmc    <- as.matrix( ModelC.coda )
  i.a            <- sapply( 1:ncv, function(i) {
    which( colnames(ModelC.mcmc)==paste("a[",i,"]",sep="") ) } )
  i.b            <- sapply( 1:ncv, function(i) {
    which( colnames(ModelC.mcmc)==paste("b[",i,"]",sep="") ) } )
  i.mu.a         <- which( colnames(ModelC.mcmc)=="mu.a"    )
  i.mu.b         <- which( colnames(ModelC.mcmc)=="mu.b"    )
  i.sigma.a      <- which( colnames(ModelC.mcmc)=="sigma.a" )
  i.sigma.b      <- which( colnames(ModelC.mcmc)=="sigma.b" )
  ModelC.a       <- ModelC.mcmc[,i.a      ] ; ModelC.b       <- ModelC.mcmc[,i.b      ]
  ModelC.mu.a    <- ModelC.mcmc[,i.mu.a   ] ; ModelC.mu.b    <- ModelC.mcmc[,i.mu.b   ]
  ModelC.sigma.a <- ModelC.mcmc[,i.sigma.a] ; ModelC.sigma.b <- ModelC.mcmc[,i.sigma.b]

# Posterior parameter distribution for the three models.
  par( mfrow=c(2,4), mar=c(3,2,3,2) )
  boxplot(ModelA.a, main="a (Model A)", ylim=c(-100,1000) )
  boxplot(ModelB.a, main="a (Model B)", names=1:6, ylim=c(-100,1000), col="cyan")
  boxplot(ModelC.a, main="a (Model C)", names=1:6, ylim=c(-100,1000), col="green")
  boxplot(list(ModelC.mu.a, ModelC.sigma.a), names=c("mu","sigma"),
          main="hyper_a (Model C)", ylim=c(-100,1000), col="green")
  boxplot(ModelA.b, main="b (Model A)", names.arg="", ylim=c(-2,4))
  boxplot(ModelB.b, main="b (Model B)", ylim=c(-2,4), col="cyan")
  boxplot(ModelC.b, main="b (Model C)", ylim=c(-2,4), col="green")
  boxplot(list(ModelC.mu.b, ModelC.sigma.b), names=c("mu","sigma"),
          main="hyper_b (Model C)", ylim=c(-2,4), col="green")
