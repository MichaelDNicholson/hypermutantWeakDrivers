#Multitype birth-death-mutation simulations
#Using W limit approximations 
#Here we assume 3 population types, precursor, with strong driver, with weak driver.


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

clone_colors = c("#619CFF",
                 "#E69F00",
                 "#F8766D" )
# 
# clone_colors = c("#619CFF",
#                  "green",
#                  "magenta" )

#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))
                 



# Clone size tracjectory plots ----------------------------------------------------------


runs <- 50
seednum <- 2
sizegrid <- 10^(seq(3,12,.2))

weak_dominatesplot <-  get_trajectory_plot(bp=0.2,
                                           dp=0,
                                           muwd=0.005,
                                           musd=0.00005,
                                           swd=1.5,
                                           ssd=2.5,
                                           seednum=seednum,
                                           runs = runs,
                                           sizegrid = sizegrid ,
                                           alp=.15)

    

strong_dominatesplot =  get_trajectory_plot(bp=0.2,
                                            dp=0,
                                            muwd=0.005/50,
                                            musd=0.00005,
                                            swd=1.5,
                                            ssd=2.5,
                                            seednum=seednum,
                                            runs = runs,
                                            sizegrid = sizegrid,
                                            alp=.15)
  
  

selectivesweep_plot =  get_trajectory_plot(bp=0.2,
                                           dp=0,
                                           muwd=0.005,
                                           musd=0.00005,
                                           swd=1.5,
                                           ssd=5,
                                           seednum=seednum,
                                           runs = runs,
                                           sizegrid = sizegrid,
                                           alp=.15)








allsim_plots = plot_grid(strong_dominatesplot,
                         weak_dominatesplot,
                selectivesweep_plot,
                labels = c('A','B','C'))

if (saveplots == T){
  plotfilename = paste0(plotdir,fileOutCommonString,
                        ".trajectory_scenarios.pdf")
  save_plot(plotfilename, allsim_plots, base_height = 4, base_asp = 1.6)
  
}


