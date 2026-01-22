
library(magrittr)


mmToInches <- function(mm){
  return(mm/25.4)
}

# relH <- mmToInches(53.7)
# relW <- mmToInches(180/3*2)
# quartz(width = relW,height = relH)



reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")

reduced_simdatadir <-paste(projectRoot,"reduced_data/simData/",sep = "")

fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOut = F


#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))






# ABC accepted selection params  -------------------------------------------


dfParamsProps <- readRDS(paste0(reduced_simdatadir,"20260108.simKrasInf50Ksims.rds"))
abcEps <- .05
dfParamsPropsAbcPass <- subset(dfParamsProps,distTotal<abcEps)




# densities of posterior params -------------------------------------------
# -----------------------------





layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
par(mar = c(4, 4, 1.2, 1),family = "Arial",cex = .6,las = 1)  

hist(
  log10(dfParamsProps$muwdNonHyp),
  breaks = "FD",
  col = rgb(0.3, 0.6, 0.9, 0.6),
  border = "white",
  xlab = expression("Non-hypermutant"~italic(mu)[weak]),
  ylab = "Counts: prior samples",
  main = "",
  yaxt = "n"
)

axis(2, las = 1)
box(bty = "l")



hist(
  log10(dfParamsPropsAbcPass$muwdNonHyp),
  breaks = "FD",
  col = rgb(0.3, 0.6, 0.9, 0.6),
  border = "white",
  xlab = expression("Non-hypermutant"~italic(mu)[weak]),
  ylab = "Counts: ABC posterior samples",
  main = "",
  yaxt = "n"
)

axis(2, las = 1)
box(bty = "l")
