# Comparison of simulations with analytic results
# Simulations use approximate form of birth-death processes

# Preamble ----------------------------------------------------------------

library(ggplot2)
library(magrittr)
library(cowplot)
library(plyr)
library(reshape2)
library(deSolve)



plotdir = paste0(projectRoot,"images/approxSimOutput/")
fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")


saveplots = F

#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))


# Contour plots: weak/strong path to third driver  -------------------------------------------------------------
# Assuming third event (driver/immune escape) occurs
# at rate proportional to weak/strong driver clone.
# Figure 4A

plotWeakGrid <- plotProb3rdEventWeakGrid(bp =.2,
                                         dp = 0,
                                         muwd = 5*10^(-3),
                                         mutBiasGrid = 10^seq(0,log10(500),by = .2),
                                         swd = 0.5,
                                         ssd = 2.5*.5,
                                         mu3grid= 10^(-(seq(8,1,-.25))),
                                         runs = 100,
                                         seednum = 2)
plotWeakGrid  <- plotWeakGrid + 
                              theme(text = element_text(family = "arial",
                                                     size = 6),
                                    legend.key.size = unit(.35,"cm"))
