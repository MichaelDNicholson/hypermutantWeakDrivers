
# Preamble ----------------------------------------------------------------

library(ggplot2)
library(magrittr)
library(cowplot)
library(plyr)
library(reshape2)
library(deSolve)


reduced_simdatadir <-paste(projectRoot,"reduced_data/simData/",sep = "")

fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOut = F


#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))


biasweakNonhyp <- 8.8
biasweakPole <- 35
mutMultPole <- 100
mu3ANonhyp <- 10^(-5)





contingencyTables <- readRDS(paste0(reduced_datadir,"contigencyTables.rds"))
contingencyTableExon2 <- contingencyTables $class1213
nPole <- contingencyTableExon2["Pole",] %>% sum
nNonHyper <- contingencyTableExon2["Non-Pole",] %>% sum
percWeakPole <- contingencyTableExon2["Pole","Weak driver"]/nPole
percWeakNonHyper<- contingencyTableExon2["Non-Pole","Weak driver"]/nNonHyper

set.seed(43)
nSamps <- 5*10^4

paramSamps <- matrix(NA,nrow = nSamps,ncol =3)
i=0
while (i<nSamps){
  muwdNonHyp <- 10^(-1*runif(1, 4,8))
  swd <- runif(1, 0,20)
  ssd <- runif(1, 0,20)
  if (swd<ssd){
    i = i+1
    paramSamps[i,] = c( muwdNonHyp, swd,ssd)
  }
}
colnames(paramSamps) = c("muwdNonHyp", "swd","ssd")

# match to observed tumour counts
prop3rdEventWeakSim = function(muwdNonhyp,swd,ssd,seedNum){
  musdNonhyp <- muwdNonhyp/biasweakNonhyp
  muwdPole <- muwdNonhyp*mutMultPole 
  musdPole <- muwdPole/biasweakPole
  
  mu3APole<- mu3ANonhyp*mutMultPole 
  probsNonHyp <- prob3rdEventWeak( bp =.2,
                                   dp = 0,
                                   muwdNonhyp ,
                                   musdNonhyp,
                                   swd,
                                   ssd,
                                   mu3ANonhyp,
                                   runs=nNonHyper,
                                   seednum=seedNum) 
  probsPole <- prob3rdEventWeak( bp =.2,
                                 dp = 0,
                                 muwdPole ,
                                 musdPole,
                                 swd,
                                 ssd,
                                 mu3APole,
                                 runs=nPole,
                                 seednum=nSamps+seedNum) 
  return(c("probsNonHyp"= probsNonHyp,
           "probsPole" = probsPole))
  
}

pb <- txtProgressBar(min = 0, max = nrow(paramSamps), style = 3)
simProps <- lapply(1:nSamps, function(k){
  setTxtProgressBar(pb, k)
  # print(k)
  prop3rdEventWeakSim( paramSamps[k,"muwdNonHyp"], 
                       paramSamps[k,"swd"],
                       paramSamps[k,"ssd"],k) %>% 
    return()
})
simProps <- do.call(rbind,simProps)

dfParamsProps <- data.frame(cbind(paramSamps , simProps))
dfParamsProps$distNonHyp = abs(dfParamsProps$probsNonHyp-percWeakNonHyper )
dfParamsProps$distPole = abs(dfParamsProps$probsPole-percWeakPole )
dfParamsProps$distTotal = dfParamsProps$distNonHyp+dfParamsProps$distPole
dfParamsProps <- dfParamsProps[order(dfParamsProps$distTotal, decreasing = F),]
dfParamsProps$delta = (dfParamsProps$ssd-dfParamsProps$swd)/((1+dfParamsProps$ssd)*(1+dfParamsProps$swd))
dfParamsProps$scaledParamNonHyp <- biasweakNonhyp/dfParamsProps$delta 
dfParamsProps$scaledParamPole <- biasweakPole/dfParamsProps$delta 

if (saveOut == T){
  filename <- paste0(reduced_simdatadir,fileOutCommonString,"simKrasInf50Ksims.rds")
  saveRDS(dfParamsProps,file = filename)
}
