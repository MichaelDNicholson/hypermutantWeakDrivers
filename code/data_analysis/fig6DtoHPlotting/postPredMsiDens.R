

library(magrittr)


mmToInches <- function(mm){
  return(mm/25.4)
}

# relH <- mmToInches(53.7)
# relW <- mmToInches(180/3)

# quartz(width = relW,height = relH)
par(mar = c(4, 4, 1.2, 1),family = "Arial",cex = .6,las = 1)  



reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")

reduced_simdatadir <-paste(projectRoot,"reduced_data/simData/",sep = "")

fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOut = F

propWeakMSIPosteriorParams <- readRDS(paste0(reduced_simdatadir,
                                             "20260109.propWeakMsiPosteriorPredSamps.rds"))

#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))


# Relevant data parameters
biasweakNonhyp <- 8.8
biasweakPole <- 35
mutMultPole <- 100
mu3ANonhyp <- 10^(-5)


dataPweakPole <- .74 #proportion non codon 12 13 in Pole tumours
dataPweakNonHyp <- .15 #proportion non codon 12 13 muts in nonhypermutant tumours




# MSI predicted proportion weak -------------------------------------------



msiWeakPropObs <- .253




dens <- density(propWeakMSIPosteriorParams)

# Plot the density
plot(dens,
     xlab =  "Proportion of weak driver MMRd CRCs " ,
     ylab = "Posterior predictive density",
     lwd = 2,        # thicker line for clarity
     col = "black",main = "")  # choose color as you like
polygon(dens, col  = rgb(0.3, 0.6, 0.9, 0.6), border = "black")
abline(v = mean(propWeakMSIPosteriorParams),lwd = 2)
abline(v = quantile(propWeakMSIPosteriorParams,.025),lwd = 2,lty = 'dashed')
abline(v = quantile(propWeakMSIPosteriorParams,.975),lwd = 2,lty = 'dashed')
abline(v= msiWeakPropObs,col = "red",lwd=2)
usr <- par("usr")  # c(xmin, xmax, ymin, ymax)
text(
  x = msiWeakPropObs,
  y = usr[4],              # top of plot
  labels = "MSI observed",
  col = "red",
  pos = 3,                 # above the point
  cex = 1,
  xpd = NA                 # allow drawing in margin if needed
)





