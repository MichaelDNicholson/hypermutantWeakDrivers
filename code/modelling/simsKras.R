
# Preamble ----------------------------------------------------------------

library(magrittr)




reduced_simdatadir <-paste(projectRoot,"reduced_data/simData/",sep = "")
fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOut = F


#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))



biasweakNonhyp <- 8.8
biasweakPole <- 35
mutMultPole <- 100


swdgrid <- seq(0.1,3,.05)
ssdgrid <- c(4,15)
runs = 5000
muwdNonhyp = 10^(-7)


pgrid <- expand.grid(swd =  swdgrid ,
                     ssd = ssdgrid)



mu3ANonhyp <- 10^(-5)
pb <- txtProgressBar(min = 0, max = nrow(pgrid), style = 3)
set.seed(43)
prob3rdEventWeakVals = sapply(1:nrow(pgrid),function(k){
  setTxtProgressBar(pb, k)
  musdNonhyp <- muwdNonhyp/biasweakNonhyp
  muwdPole <- muwdNonhyp*mutMultPole 
  musdPole <- muwdPole/biasweakPole
 
  mu3APole<- mu3ANonhyp*mutMultPole 
  probsNonHyp <- prob3rdEventWeak( bp =.2,
                    dp = 0,
                    muwdNonhyp ,
                    musdNonhyp,
                    pgrid[k,"swd"],
                    pgrid[k,"ssd"],
                    mu3ANonhyp,
                    runs=runs ,
                    seednum=k) 
  probsPole <- prob3rdEventWeak( bp =.2,
                                   dp = 0,
                                   muwdPole ,
                                   musdPole,
                                   pgrid[k,"swd"],
                                   pgrid[k,"ssd"],
                                   mu3APole,
                                   runs=runs ,
                                   seednum=nrow(pgrid)+k) 
  return(c("probsNonHyp"= probsNonHyp,
           "probsPole" = probsPole))

}) 

dfProb3rdPgrid <- data.frame(cbind(t(prob3rdEventWeakVals),pgrid))



if (saveOut == T){
  filename <- paste0(reduced_simdatadir ,fileOutCommonString,"simOutWeakPropVsWeakSel.rds")
  saveRDS(dfProb3rdPgrid,file = filename)
}



# Prop vs weak bias Kras simulations ----------------------------------------------------


dfParamsProps <- readRDS(paste0(reduced_simdatadir,"20260108.simKrasInf50Ksims.rds"))
abcEps <- .05
dfParamsPropsAbcPass <- subset(dfParamsProps,distTotal<abcEps)


muwdNonhyp  <- dfParamsPropsAbcPass[1,"muwdNonHyp"] 
swd  <- dfParamsPropsAbcPass[1,"swd"]
ssd  <- dfParamsPropsAbcPass[1,"ssd"]
biasweakNonhyp <- 8.8
mu3ANonhyp <- 10^(-5)
biasweakPoleGrid <- 10^seq(0,2,.2)
mutMultPoleGrid <- c(1,5,20,50,100)

pgrid = expand.grid(biasweakPole = biasweakPoleGrid ,
                    mutMultPole=  mutMultPoleGrid)
mu3ANonhyp <- 10^(-5)
pb <- txtProgressBar(min = 0, max = nrow(pgrid), style = 3)
nRuns <- 50000
set.seed(43)
prob3rdEventWeakVals = sapply(1:nrow(pgrid),function(k){
  setTxtProgressBar(pb, k)
  
  musdNonhyp <- muwdNonhyp/biasweakNonhyp
  
  muwdPole <- muwdNonhyp*pgrid[k,"mutMultPole"]
  musdPole <- muwdPole/pgrid[k,"biasweakPole"]
  
  mu3APole<- mu3ANonhyp*pgrid[k,"mutMultPole"]
  
  probsNonHyp <- prob3rdEventWeak( bp =.2,
                                   dp = 0,
                                   muwdNonhyp ,
                                   musdNonhyp,
                                   swd,
                                   ssd,
                                   mu3ANonhyp,
                                   runs=nRuns,
                                   seednum=k) 
  probsPole <- prob3rdEventWeak( bp =.2,
                                 dp = 0,
                                 muwdPole ,
                                 musdPole,
                                 swd,
                                 ssd,
                                 mu3APole,
                                 runs=nRuns,
                                 seednum=nrow(pgrid)+k) 
  return(c("probsNonHyp"= probsNonHyp,
           "probsPole" = probsPole))
  
}) 



dfProb3rdPgridWeakBias <- data.frame(cbind(t(prob3rdEventWeakVals),pgrid))


if (saveOut == T){
  filename <- paste0(reduced_simdatadir ,fileOutCommonString,"simOutWeakPropVsWeakbias.rds")
  saveRDS(dfProb3rdPgridWeakBias,file = filename)
}


# MSI predicted proportion weak -------------------------------------------


biasweakMSI <- 6.73
mutMultMSI <- 10
nMSI = 82+68 #number of MSI in cbio + number of MSI in 100kgp
msiWeakPropObs <- .25

set.seed(43)



# match to observed tumour counts
    prop3rdEventWeakMSISim = function(muwdNonhyp,swd,ssd,mutMultMSI,seedNum){
      muwdMSI <- muwdNonhyp*mutMultMSI
      musdMSI <- muwdMSI/biasweakMSI
      mu3ANonhyp <- 10^(-5)
      mu3AMSI<- mu3ANonhyp*mutMultMSI
      probsMSI<- prob3rdEventWeak( bp =.2,
                                   dp = 0,
                                   muwdMSI ,
                                   musdMSI ,
                                   swd,
                                   ssd,
                                   mu3AMSI,
                                   runs=nMSI,
                                   seednum=seedNum) 
      
      return(c("probsMSI"= probsMSI))
      
    } 

# propWeakByMsiMult <- lapply(seq(5,50,5), function(mutMultMSI){

    


  propWeakMSIPosteriorParams <- sapply(1:nrow(dfParamsPropsAbcPass), function(k){
    prop3rdEventWeakMSISim (dfParamsPropsAbcPass$muwdNonHyp[k],
                            dfParamsPropsAbcPass$swd[k],
                            dfParamsPropsAbcPass$ssd[k],
                            mutMultMSI,
                            k)
  })



if (saveOut == T){
  filename <- paste0(reduced_simdatadir ,fileOutCommonString,"propWeakMsiPosteriorPredSamps.rds")
  saveRDS(propWeakMSIPosteriorParams,file = filename)
}



   
