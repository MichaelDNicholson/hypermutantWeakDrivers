

# Theory functions ---------------------------------------------------------



#### Theory: switch subtypes ####
theoryRegimeFun <- function(biasA,biasB,mu3A,mu3B,del,
                            general=T){
  
  #biasA = muwdA/musdA
  # biasB =  muwdB/musdB
  # mu3A = mu3A
  # mu3B = mu3B
  # del =   delGrid[1]
  
  x = mu3B/mu3A
  weakA <- (biasA>(mu3A^(-del))) #weak seen in subtype A
  weakB <- (biasB>(mu3B^(-del))) #weak seen in subtype B
  
  
  if(general == T){
    
    if (weakA & weakB){
      return("Subtype A: Weak driver\nSubtype B: Weak driver")
    } else if ( !weakA & weakB){
      return("Subtype A: Strong driver\nSubtype B: Weak driver")
      # } else if (weakA & !weakB){
      #assume muweak>mustrong
      #   return("Subtype A: Weak driver\nPOLE CRC: Strong driver")
    } else{
      return("Subtype A: Strong driver\nSubtype B: Strong driver")
    }
    
  } else {
    
    if (weakA & weakB){
      return("Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Weak KRAS")
    } else if ( !weakA & weakB){
      return("Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Weak KRAS")
      # } else if (weakA & !weakB){
      #   return("Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Strong KRAS")
      #assume muweak>mustrong
    } else{
      return("Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Strong KRAS")
    }
    
  }
}

mainRegimes <- c("Subtype A: Strong driver\nSubtype B: Strong driver",
                 "Subtype A: Strong driver\nSubtype B: Weak driver",
                 "Subtype A: Weak driver\nSubtype B: Weak driver")

mainRegimesSpecific <- c("Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Strong KRAS",
                         "Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Weak KRAS",
                         "Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Weak KRAS")
# colfunc <- colorRampPalette(c("lightgoldenrod", "red4"))
mainRegimeColors <- c("red4","darkorange","lightgoldenrod1") 

biasChangeString <- expression(paste("Fold increase in weak driver bias (",beta["B"],'/',
                                     beta["A"],')'))
mutChangeString <- expression(paste("Fold increase in tertiary driver rate (",mu["3"]^"B",'/',
                                    mu["3"]^"A",')'))

biasChangeString <- expression(paste("Fold increase in weak driver bias (",beta^"B",'/',
                                     beta^"A",')'))

# deltaSelString <- expression(atop('Scaled selection difference',
#                                   paste(Delta,'=(',s["strong"]-s["weak"],')/(1+',s["strong"],
#                                         ')(1+',s["weak"],
#                                         ')')))
deltaSelString <- expression(paste("Scaled selection difference (",Delta,')'))


# General analytic regimes plots ------------------------------------------
# Figure 4C


plotContourTheory_mutMultGridBiasMultGrid_SwitchThirdRouteContour <- function(biasA,
                                                                              biasMultGrid,
                                                                              del,
                                                                              mu3A,
                                                                              mutMultGrid,
                                                                              general){
  
  # 
  # del = .2
  # mu3A = .000001
  # biasA <- 2
  # 
  # general = T
  # # # mutMultGrid <- seq(10,250,5)
  # mutMultGrid <- 10^(seq(0,log10(500),.2))
  # biasMultGrid <-  10^(seq(0,log10(500),.2))
  # 
  
  pgrid <- expand.grid("bmult" =  biasMultGrid ,
                       "x" =  mutMultGrid)
  
  
  
  
  
  
  
  theoryRegime <- sapply(1:nrow(pgrid), function(k){
    theoryRegimeFun(biasA,
                    biasA*pgrid[k,"bmult"],
                    mu3A,
                    pgrid[k,"x"]*mu3A,
                    del) %>% return()
  })
  
  # mainRegimes <- c("Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Strong KRAS",
  #                  "Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Weak KRAS",
  #                  "Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Weak KRAS")
  
  regimeTheory <- factor(theoryRegime,
                         levels = c(mainRegimes))
  
  switch3rdEvent<- cbind(pgrid,
                         "regimeTheory" =  regimeTheory) %>% 
    as.data.frame()
  
  
  
  
  switch3rdEvent$z = rep(1,nrow(switch3rdEvent))
  switch3rdEvent$z[switch3rdEvent$regimeTheory == mainRegimes[1]] =1                             
  switch3rdEvent$z[switch3rdEvent$regimeTheory == mainRegimes[2]] =2                             
  switch3rdEvent$z[switch3rdEvent$regimeTheory == mainRegimes[3]] =3                             
  
  
  pltcontourTheory <-  ggplot()+
    geom_contour_filled(data =   switch3rdEvent,
                        aes(x= x,y=  bmult,
                            z=z),
                        breaks = seq(.5,3.5,1))+
    scale_y_continuous(trans = "log10")+ 
    scale_x_continuous(trans = "log10" )+
    ylab(   biasChangeString )+
    xlab(  mutChangeString)+
    theme_bw()+
    theme(legend.position = "right",
          legend.spacing.y = unit(.1, 'cm'))+
    guides(fill = guide_legend(byrow = TRUE,
                               ncol = 1))+
    scale_fill_manual(values= mainRegimeColors,
                      name = "Expected driver class\nacross tumour subtypes",
                      labels = mainRegimes,
                      drop=FALSE)
  
  # pltcontourTheory 
  
  
  
  return(  pltcontourTheory )
}

plotContourTheory_delMu3MultGrid_SwitchThirdRouteContour <- function(biasA,
                                                                     biasB,
                                                                     swd,
                                                                     delGridDelta,
                                                                     mu3A,
                                                                     mutMultGrid,
                                                                     general){
  
  # biasA = 22
  # biasB = 44
  # swd=.3
  # delGridDelta= .01
  # mu3A = .00001
  # mutMultGrid= 10^(seq(0,2.5,.01))
  # general = T
  # 
  # 
  # 
  # #if swd = 0.35, then -delta =  (1/(1+swd)-1/(1+ssd))
  # # = 0.74-1/(1+ssd)-> max is 0.74
  # #(weirdly stronger swd then less values possible for #ssd
  # # such that condition is true
  # swd = .3
  # delGridDelta = .01
  # mu3A = .000001
  # 
  # 
  # # mutMultGrid <- seq(10,250,5)
  # mutMultGrid <- 10^(seq(0,log10(500),.2))
  
  delGrid <- seq(delGridDelta,1/(1+swd),delGridDelta)
  
  pgrid <- expand.grid("del" =  delGrid,
                       "x" =  mutMultGrid)
  
  
  
  
  
  
  theoryRegime <- sapply(1:nrow(pgrid), function(k){
    theoryRegimeFun(biasA,
                    biasB,
                    mu3A,
                    pgrid[k,"x"]*mu3A,
                    pgrid[k,"del"]) %>% return()
  })
  
  # mainRegimes <- c("Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Strong KRAS",
  #                  "Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Weak KRAS",
  #                  "Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Weak KRAS")
  
  regimeTheory <- factor(theoryRegime,
                         levels = c(mainRegimes))
  
  switch3rdEventSelection <- cbind(pgrid,
                                   "regimeTheory" =  regimeTheory) %>% 
    as.data.frame()
  
  switch3rdEventSelection$z = rep(1,nrow(switch3rdEventSelection))
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == mainRegimes[1]] =1                             
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == mainRegimes[2]] =2                             
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == mainRegimes[3]] =3                             
  
  
  
  
  # atop("Mutational bias to weak driver",
  #      paste('(', mu["weak"]/mu["strong"], ')')))
  # 
  # limits = c(min(switch3rdEventSelection$x),
  #            max(switch3rdEventSelection$x))
  
  pltcontourTheory <-  ggplot()+
    geom_contour_filled(data =   switch3rdEventSelection,
                        aes(x= x,y= del,
                            z=z),
                        breaks = seq(.5,3.5,1))+
    scale_y_continuous(trans = "identity")+ 
    scale_x_continuous(trans = "log10" )+
    scale_fill_manual(values= mainRegimeColors,
                      name = "Expected driver class\nacross tumour subtypes",
                      labels = mainRegimes,
                      drop=FALSE)+
    ylab(  deltaSelString )+
    xlab(  mutChangeString)+
    theme_bw()+
    theme(legend.position = "right",
          legend.spacing.y = unit(.1, 'cm'))+
    guides(fill = guide_legend(byrow = TRUE,
                               ncol = 1))
  
  # pltcontourTheory 
  
  
  
  return(  pltcontourTheory )
}

plotContourTheory_biasdel_SwitchThirdRouteContour <- function(biasA ,
                                                              biasMultGrid, 
                                                              swd,
                                                              delGridDelta,
                                                              mu3A,
                                                              mutMult,
                                                              general){
  
  
  # biasA <- 10
  # biasMultGrid <-  10^(seq(0,log10(500),.02))
  # swd = .3
  # delGridDelta <-  .01
  # mu3A <-  .000001
  # mutMult <- 10
  # 
  
  
  delGrid <- seq(delGridDelta,1/(1+swd),delGridDelta)
  mu3B <- mutMult*mu3A
  
  pgrid <- expand.grid("del" =  delGrid,
                       "bmult" =  biasMultGrid )
  
  
  
  
  
  
  theoryRegime <- sapply(1:nrow(pgrid), function(k){
    theoryRegimeFun(biasA,
                    biasA*pgrid[k,"bmult"],
                    mu3A,
                    mu3B,
                    pgrid[k,"del"]) %>% return()
  })
  
  # mainRegimes <- c("Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Strong KRAS",
  #                  "Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Weak KRAS",
  #                  "Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Weak KRAS")
  
  regimeTheory <- factor(theoryRegime,
                         levels = c(mainRegimes))
  
  switch3rdEventSelection <- cbind(pgrid,
                                   "regimeTheory" =  regimeTheory) %>% 
    as.data.frame()
  
  switch3rdEventSelection$z = rep(1,nrow(switch3rdEventSelection))
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == mainRegimes[1]] =1                             
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == mainRegimes[2]] =2                             
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == mainRegimes[3]] =3                             
  
  
  
  
  
  pltcontourTheory <-  ggplot()+
    geom_contour_filled(data =   switch3rdEventSelection,
                        aes(x= bmult,y= del,
                            z=z),
                        breaks = seq(.5,3.5,1))+
    scale_y_continuous(trans = "identity")+ 
    scale_x_continuous(trans = "log10" )+
    scale_fill_manual(values= mainRegimeColors,
                      name = "Expected driver class\nacross tumour subtypes",
                      labels = mainRegimes,
                      drop=FALSE)+
    ylab(  deltaSelString )+
    xlab(  biasChangeString )+
    theme_bw()+
    theme(legend.position = "right",
          legend.spacing.y = unit(.1, 'cm'))+
    guides(fill = guide_legend(byrow = TRUE,
                               ncol = 1))
  # ggtitle(mu3String )
  
  # pltcontourTheory 
  
  
  
  return(  pltcontourTheory )
}


# KRAS POLE regimes -------------------------------------------------------
# Fig 6


regimeLabels <- c("Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Weak KRAS",
                  "Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Weak KRAS",
                  "Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Strong KRAS")

plotProb3rdEventCancer_delMut3AGridContour <-   function(biasA,
                                                         biasB,
                                                         swd,
                                                         delGridDelta,
                                                         mu3Agrid,
                                                         mutMult,
                                                         general){
  
  
  # 
  # 
  # 
  # #if swd = 0.35, then -delta =  (1/(1+swd)-1/(1+ssd))
  # # = 0.74-1/(1+ssd)-> max is 0.74
  # #(weirdly stronger swd then less values possible for #ssd
  # # such that condition is true
  # swd = .3
  # delGridDelta = .01
  # biasA <- 20
  # biasB <- 40
  # #
  # #
  # mutMult <- 100
  # mu3Agrid <- 10^seq(-9,-3,.5)
  # general = F
  
  delGrid <- seq(delGridDelta,1/(1+swd),delGridDelta)
  
  pgrid <- expand.grid("del" =  delGrid,
                       "mu3A" =   mu3Agrid)
  
  
  
  
  
  
  theoryRegime <- sapply(1:nrow(pgrid), function(k){
    theoryRegimeFun(biasA,
                    biasB,
                    pgrid[k,"mu3A"],
                    pgrid[k,"mu3A"]*mutMult,
                    pgrid[k,"del"], 
                    general) %>% return()
  })
  
  # mainRegimes <- c("Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Strong KRAS",
  #                  "Nonhypermutant CRC: Strong KRAS\nPOLE CRC: Weak KRAS",
  #                  "Nonhypermutant CRC: Weak KRAS\nPOLE CRC: Weak KRAS")
  if (general == T){
    regimeTheory <- factor(theoryRegime,
                           levels = c(mainRegimes))
  } else {
    regimeTheory <- factor(theoryRegime,
                           levels = c(mainRegimesSpecific))
  }
  
  switch3rdEventSelection <- cbind(pgrid,
                                   "regimeTheory" =  regimeTheory) %>% 
    as.data.frame()
  
  
  
  
  mutAString <- expression(paste("Tertiary driver rate in\nnon-hypermutant CRC (",
                                 mu["3"]^"A",')'))
  
  switch3rdEventSelection$z = rep(1,nrow(switch3rdEventSelection))
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == regimeLabels[1]] =1                             
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == regimeLabels[2]] =2                             
  switch3rdEventSelection$z[switch3rdEventSelection$regimeTheory == regimeLabels[3]] =3                             
  
  
  pltcontourTheory <-  ggplot()+
    geom_contour_filled(data =   switch3rdEventSelection,
                        aes(x= mu3A,y= del,
                            z=z),
                        breaks = seq(.5,3.5,1))+
    scale_y_continuous(trans = "identity")+ 
    scale_x_continuous(trans = "log10" )+
    scale_fill_manual(values= mainRegimeColors[3:1],
                      label = regimeLabels,
                      name = "Expected driver class\nacross tumour subtypes",
                      drop=FALSE)+
    ylab(  deltaSelString )+
    xlab(  mutAString)+
    theme_bw()+
    theme(legend.position = "right",
          legend.spacing.y = unit(.1, 'cm'))+
    guides(fill = guide_legend(byrow = TRUE,
                               ncol = 1))
  
  # pltcontourTheory 
  
  
  
  return(  pltcontourTheory )
}


plotSelRegimeDriverFixedMut3AContour <- function(biasweakNonhyp,
                                                 biasweakPole,
                                                 mu3A,
                                                 mutMultPole,
                                                 swdgrid,
                                                 ssdgrid){
  
  pgrid <- expand.grid("swd" = swdgrid,
                       "ssd" = ssdgrid)
  theoryRegime <- sapply(1:nrow(pgrid), function(k){
    
    swd = pgrid[k,"swd"]
    ssd = pgrid[k,"ssd"]
    if (swd>ssd){
      return("NA")
    } else {
      del <- delta_fun(swd,ssd)
      theoryRegimeFun(biasweakNonhyp,
                      biasweakPole,
                      mu3A,
                      mu3A*mutMultPole,
                      del,
                      F) %>% return()
    }
  })
  #NA because swd has to be <= ssd
  regimeTheory <- factor(theoryRegime,
                         levels = c(mainRegimesSpecific,"NA"))
  
  regimeTheorydf <- cbind(pgrid,regimeTheory) %>% as.data.frame()
  regimeTheorydOuts <- regimeTheorydf$regimeTheory %>% unique
  regimeTheorydf$z <- rep(0,nrow(regimeTheorydf))
  regimeTheorydf$z[regimeTheorydf$regimeTheory=="NA"] = 0
  regimeTheorydf$z[regimeTheorydf$regimeTheory== regimeTheorydOuts[1]] = 1
  regimeTheorydf$z[regimeTheorydf$regimeTheory==regimeTheorydOuts[3]] = 2
  regimeTheorydf$z[regimeTheorydf$regimeTheory==regimeTheorydOuts[4]] = 3
  
  
  
  pltSvalForFixedU3A <- ggplot()+
    geom_contour_filled(data =   subset(regimeTheorydf,
                                        swd<=ssd),
                        aes(x= swd,y= ssd,
                            z=z),
                        breaks = seq(.5,3.5,1))+
    xlim(c(0,3))+
    scale_fill_manual(values= c(mainRegimeColors[3:1]),
                      name = "Expected driver class\nacross tumour subtypes",
                      labels=regimeLabels ,
                      drop=FALSE)+
    theme_bw()+
    xlab(expression(paste("Weak KRAS selection coefficient (",s["weak"],')')))+
    ylab(expression(paste("Strong KRAS selection coefficient (",s["strong"],')')))+
    theme(legend.position = "none")
  
  return(pltSvalForFixedU3A)
}
