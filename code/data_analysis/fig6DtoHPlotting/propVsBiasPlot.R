

library(magrittr)


mmToInches <- function(mm){
  return(mm/25.4)
}

# relH <- mmToInches(53.7)
# relW <- mmToInches(180/3)
# quartz(width = relW,height = relH)

par(mar = c(4, 4, 1, 1),family = "Arial",cex = .6,las = 1)  

reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")

reduced_simdatadir <-paste(projectRoot,"reduced_data/simData/",sep = "")

fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOut = F


#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))


# Relevant data parameters
biasweakNonhyp <- 8.8
biasweakPole <- 35
mutMultPole <- 100
mu3ANonhyp <- 10^(-5)


dataPweakPole <- .74 #proportion non codon 12 13 in Pole tumours
dataPweakNonHyp <- .15 #proportion non codon 12 13 muts in nonhypermutant tumours

dfProb3rdPgridWeakBias<- readRDS(paste0(reduced_simdatadir,
                                    "20260114.simOutWeakPropVsWeakbias.rds"))


# Proportion of weak KRAS vs mu bias  -------------------------------------


dataPweakPole <- .74 #proportion non codon 12 13 in Pole tumours
dataPweakNonHyp <- .15 #proportion non codon 12 13 muts in nonhypermutant tumours

# What is the below plot? 
# fix a muwd in nonHyp, then all rest is fixed, then fix strong driver select coef
# and find where non-hyp and 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7","purple")

plot(subset(dfProb3rdPgridWeakBias,mutMultPole==1)$biasweakPole,
     subset(dfProb3rdPgridWeakBias,mutMultPole==1)$probsPole,ylim = c(0,.85),xlim = c(1,100),
     type = 'l',col =cbbPalette[2],lwd = 2,
     xlab = expression(paste(italic("POLE"), "-mutant (", mu[weak] / mu[strong], ")")),
     ylab = "Tumour proportion with weak driver",log = 'x')
# mtext( "Tumour proportion with weak driver",
#        side = 2,
#        line = 2,
#        adj = 1,
#        cex = .6) 
lines(subset(dfProb3rdPgridWeakBias,mutMultPole==5)$biasweakPole,
      subset(dfProb3rdPgridWeakBias,mutMultPole==5)$probsPole,
      col = cbbPalette[3],lwd = 2)

lines(subset(dfProb3rdPgridWeakBias,mutMultPole==20)$biasweakPole,
      subset(dfProb3rdPgridWeakBias,mutMultPole==20)$probsPole,
      col = cbbPalette[4],lwd = 2)


lines(subset(dfProb3rdPgridWeakBias,mutMultPole==100)$biasweakPole,
      subset(dfProb3rdPgridWeakBias,mutMultPole==100)$probsPole,
      col = cbbPalette[8],lwd = 2)

legend(x = .7, y = 0.9,
       legend = c("1X", "5X", "20X", "100X"),
       col = c(cbbPalette[2], cbbPalette[3], cbbPalette[4], cbbPalette[8]),
       lwd = 2,
       bty = "n",
       text.width = strwidth("100X"),  # ensures alignment
       title = "Fold increase muts/div",
       title.adj = 2,
       x.intersp = 0.5)
