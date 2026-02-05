# Plots exploring analytic regimes

# Preamble ----------------------------------------------------------------

library(ggplot2)
library(magrittr)
library(cowplot)
library(plyr)
library(reshape2)



plotdir = paste0(projectRoot,"images/paramRegimesAnalytic/")
reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")

saveplots = F
fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")


source(paste0(projectRoot,"code/modelling/theoryFunctionsClean.R"))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                         




# Delta function vs sel parameters ----------------------------------------
# Fig 4B

delta_fun <- function(swd,ssd){
  
  return(( ssd-swd)/((1+swd)*(1+ssd)))
}
swdgrid <- c(0.2,0.5,1,2)
ssdgrid <- seq(.21,20,.1)

delGridValsSwd <- sapply(swdgrid, function(swd){
  sapply(ssdgrid , function(x) delta_fun(swd,x)) %>% return()
})

delGridValsSwdSsdGrid <- cbind(ssdgrid,delGridValsSwd) %>% as.data.frame()


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                         

delGridValsSwdSsdGrid_long <- melt(delGridValsSwdSsdGrid,
                                   id = c("ssdgrid"))
# deltaSelString <- expression(atop('Scaled selection difference',
#                                   paste(Delta,'=(',s["strong"]-s["weak"],')/(1+',s["strong"],
#                                         ')(1+',s["weak"],
#                                         ')')))
deltaSelString <- expression("Scaled selection difference"~(Delta))

plotDelFn <- ggplot(delGridValsSwdSsdGrid_long,aes(x=ssdgrid,
                                                   y = value,
                                                   color = variable))+
  geom_line(size = .5)+
  theme_bw(base_family = 'arial')+
  ylim(c(0,max(max(delGridValsSwd))))+
  scale_color_manual(values = cbbPalette[1:4],
                     name= expression(s["weak"]),
                     labels = as.character(swdgrid))+
  xlab(expression(s["strong"]))+
  ylab(deltaSelString )+
  theme(text = element_text(family = "arial",size = 6))


if (saveplots == T){
  pltfilename = paste0(plotdir, fileOutCommonString,
                       "plotDelFn.pdf")
 
  ggplot2::ggsave(filename =   pltfilename , 
                  plot =plotDelFn  , 
                  device = cairo_pdf, 
                  dpi = 200, 
                  width = 6,
                  height = 4, 
                  units = "cm")
  
  
  
  
}



# Switch driver between subtypes ------------------------------------------
# Fig 4C


#common params
biasA <- 10
fixdel <- 0.5
mu3A <- 10^(-5)
biasB <- 20
swd <- .3
mutMult <- 10

plotMutModBiasModSubtypes <- plotContourTheory_mutMultGridBiasMultGrid_SwitchThirdRouteContour(biasA=biasA,
                                                           biasMultGrid = 10^(seq(0,log10(100),.02)),
                                                           del = fixdel,
                                                           mu3A = mu3A,
                                                           mutMultGrid = 10^(seq(0,log10(100),.02)),
                                                           general = T)+
                              theme(legend.position = "none")+
                            theme(text = element_text(family = "arial",size = 6))

  


plotDelMutModSubtypes <- plotContourTheory_delMu3MultGrid_SwitchThirdRouteContour(biasA = biasA,
                                                  biasB = biasB,
                                                  swd = swd,
                                                  delGridDelta= .01,
                                                  mu3A = mu3A,
                                                  mutMultGrid= 10^(seq(0,2.5,.01)),
                                                  general = T)+
                                                  theme(legend.position = "none")+
                                                    theme(text = element_text(family = "arial",size = 6))


plotBiasModDelSubtypes <- plotContourTheory_biasdel_SwitchThirdRouteContour(biasA =biasA,
                                                       biasMultGrid=10^(seq(0,log10(100),.02)), 
                                                       swd= swd,
                                                       delGridDelta=.01,
                                                       mu3A = mu3A,
                                                       mutMult <- mutMult ,
                                                       general = T)+
                                                        theme(text = element_text(family = "arial",size = 6))

relLeg <- get_legend(plotBiasModDelSubtypes)
plotBiasModDelSubtypes <- plotBiasModDelSubtypes+ theme(legend.position = "none")
allPlotSwitchDriver <- plot_grid(plotMutModBiasModSubtypes ,
          plotDelMutModSubtypes,
          plotBiasModDelSubtypes,
          relLeg,
          nrow = 1,
          rel_widths = c(1,1,1,.75))





# biasA <- 10
# biasMultGrid <-  10^(seq(0,log10(500),.02))
# swd = .3
# delGridDelta <-  .01
# mu3A <-  .000001
# mutMult <- 10


if (saveplots == T){
  

  pltfilename = paste0(plotdir,fileOutCommonString,
                       "allSwitchDriver.pdf")
  ggplot2::ggsave(filename =   pltfilename , 
                  plot =  allPlotSwitchDriver , 
                  device = cairo_pdf, 
                  dpi = 250, 
                  width = 18,
                  height = 4.5, 
                  units = "cm")

  

  
}


# POLE KRAS parameter regimes ---------------------------------------------
# Fig 6 A-C

possibleClassVarsExplanation <- c(
  "Bioactivation score",
  "Exon 2 mutations are strong drivers",
  "observed:expected",
  "observed:expected, MSS only"
  
)



getExpectedDriverPolePlots <- function(biasweakNonhyp,
                                       biasweakPole,
                                       mutMultPole,
                                       swd,
                                       mu3A,
                                       swdgrid,
                                       ssdgrid){


  
 
    plotDelMut3APOLEvals <-  plotProb3rdEventCancer_delMut3AGridContour(biasA = biasweakNonhyp,
                                                                 biasB = biasweakPole,
                                                                 swd = swd ,
                                                                 delGridDelta = .01,
                                                                 mu3Agrid = 10^seq(-9,-3,.01),
                                                                 mutMult = mutMultPole,
                                                                 general =F)+
      theme(text = element_text(family = "arial",size = 6),
            axis.title.x = element_text(vjust=-8))
    
    plotDelMut3APOLEvalsnoMult <-  plotProb3rdEventCancer_delMut3AGridContour(biasA = biasweakNonhyp,
                                                                       biasB =biasweakPole,
                                                                       swd = swd,
                                                                       delGridDelta = .01,
                                                                       mu3Agrid = 10^seq(-9,-3,.02),
                                                                       mutMult = 1,
                                                                       general =F)+
      theme(text = element_text(family = "arial",size = 6),
            legend.position = "none",
           axis.title.x = element_text(vjust=-8))
    
    
    
    pltSvalForFixedU3A <- plotSelRegimeDriverFixedMut3AContour(biasweakNonhyp,
                                                        biasweakPole,
                                                        mu3A,
                                                        mutMultPole,
                                                        swdgrid,
                                                        ssdgrid) +
      theme(text = element_text(family = "arial",size = 6))
    
    relLeg <- get_legend(plotDelMut3APOLEvals)
    
    plotDelMut3APOLEvals <- plotDelMut3APOLEvals+ theme(legend.position = "none")
    allPlotSwitchDriverPole <- plot_grid(plotDelMut3APOLEvals ,
                                     plotDelMut3APOLEvalsnoMult ,
                                     pltSvalForFixedU3A ,
                                     relLeg,
                                     nrow = 1,
                                     rel_widths = c(1,1,1,.75))


    
  return(allPlotSwitchDriverPole)
}

inferredBiases <- read.csv(file = paste0(reduced_datadir,"ratioMutBiasPoleNonHyp.csv"))

# codon 12/13 params:
biasweakNonhyp <- 8.8
biasweakPole <- 35
mutMultPole <- 100
swd <- .3

mu3A <- 10^(-5) 
swdgrid <- seq(0.1,3,.02)
ssdgrid <- seq(0.1,20,.025)

allPlotSwitchDriverPoleCriteria <- lapply(1:nrow(inferredBiases), function(k){
  
  allPlotSwitchDriverPole <-getExpectedDriverPolePlots(inferredBiases$WeakBiasNonHyp[k],
                                                       inferredBiases$WeakBiasPole[k],
                                                       mutMultPole,
                                                       swd,
                                                       mu3A,
                                                       swdgrid,
                                                       ssdgrid)
  
  classTitle <- ggdraw() + 
    draw_label(possibleClassVarsExplanation[k],
               fontface='bold')
  allPlotSwitchDriverPoleTitle <- plot_grid(classTitle,
                                            allPlotSwitchDriverPole, nrow =2,
                                            rel_heights = c(.2,1))
  return(allPlotSwitchDriverPoleTitle)
  
})

allPlotSwitchDriverPoleCodon1213 <- allPlotSwitchDriverPoleCriteria[[2]]

allPlotSwitchDriverPoleAllCrit <- plot_grid(plotlist = allPlotSwitchDriverPoleCriteria,ncol=1)
if (saveplots == T){
  
  mmToInch <- function(m){return(m/25.4)}
  
  pltfilename = paste0(plotdir,fileOutCommonString,
                       "allSwitchDriverPoleContour.pdf")
  ggplot2::ggsave(filename =   pltfilename , 
                  plot =  allPlotSwitchDriverPoleCodon1213, 
                  device = cairo_pdf, 
                  dpi = 250, 
                  width = 18,
                  height = 4.5, 
                  units = "cm")
  
  pltfilename = paste0(plotdir,fileOutCommonString,
                       "allSwitchDriverPoleAllCrit.pdf")
  ggplot2::ggsave(filename =   pltfilename , 
                  plot =  allPlotSwitchDriverPoleAllCrit, 
                  device = cairo_pdf, 
                  dpi = 250, 
                  width = 18,
                  height = 4.5*4, 
                  units = "cm")
  
  
}
