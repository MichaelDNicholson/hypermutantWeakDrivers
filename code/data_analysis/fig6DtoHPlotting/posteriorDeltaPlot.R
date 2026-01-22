

library(magrittr)


mmToInches <- function(mm){
  return(mm/25.4)
}

# relH <- mmToInches(53.7)
# relW <- mmToInches(180/3)
# 
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

# ABC accepted selection params  -------------------------------------------


dfParamsProps <- readRDS(paste0(reduced_simdatadir,"20260108.simKrasInf50Ksims.rds"))
abcEps <- .05
dfParamsPropsAbcPass <- subset(dfParamsProps,distTotal<abcEps)





# -----------------------------

# delta




dens <- density(dfParamsPropsAbcPass$delta)

# Plot the density
plot(dens,
     xlab =   expression("Scaled selection difference (" * Delta * ")"),
     ylab = "Posterior density",
     lwd = 2,        # thicker line for clarity
     col = "black",main = "")  # choose color as you like
polygon(dens,   col = rgb(0.3, 0.6, 0.9, 0.6), border = "black")
abline(v = mean(dfParamsPropsAbcPass$delta),lwd = 2)
abline(v = quantile(dfParamsPropsAbcPass$delta,.025),lwd = 2,lty = 'dashed')
abline(v = quantile(dfParamsPropsAbcPass$delta,.975),lwd = 2,lty = 'dashed')

