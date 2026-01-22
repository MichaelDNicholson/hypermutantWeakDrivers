# Script to annotate Kras mutations as weak/strong drivers
# Will add columns with weak/strong under varied criteria.
    
# Preamble ----------------------------------------------------------------
library(readxl)
library(plyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)


  "/Users/mnichol3/Library/CloudStorage/OneDrive-UniversityofEdinburgh/oneDriveMN/",
  "research_projects/XDF_rotation_projects/weakDriversClean/"
)
reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")
miscelDir <- paste(projectRoot,"data/miscellaneous/",sep = "")
plotdir = paste0(projectRoot,"images/krasPublicData/")
fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOutput <- F

#for plotting
clone_colors = c("#619CFF",
                "darkgoldenrod1",
                "#F8766D" )

#load in dataframe of kras mutations
combinedKrasList <- read.csv(file =paste0(reduced_datadir,"combinedKrasList.csv"),
                             header = T)




# Haigis Bioactivation score ----------------------------------------------


# Use scores from Fig 3E  https://doi.org/10.1016/j.trecan.2017.08.006
# Mutations without activation scores classified as weak.s
haigisActivationScore <- read_xlsx(path = paste0(miscelDir, 
                                                 "Haigis18Activationscore.xlsx"),
                         skip =1)
haigisActivationScore$mutAA <- sapply(haigisActivationScore$...1,function(x){
  paste0("p.",x) %>% return()
})

combinedKrasList$haigisActivationScore <- sapply(combinedKrasList$mutAA, function(x){
  if (x %in% haigisActivationScore$mutAA){
    return(subset(haigisActivationScore, mutAA == x)$`Haigis Activation Score`)
  } else {
    return(NA)
  }
})

actScores <- na.omit(combinedKrasList$haigisActivationScore)
kmeanSplitActScore <- actScores %>% 
                      kmeans(.,2)
actClass <- rep("Weak",length(actScores))
if (which.max(kmeanSplitActScore$centers)==2){
  #2nd cluster is biggest
  actClass[kmeanSplitActScore$cluster==2] = "Strong"
} else {
  actClass[kmeanSplitActScore$cluster==1] = "Strong"
}

combinedKrasList$classActScore <- rep("Weak", nrow(combinedKrasList))
combinedKrasList$classActScore [ combinedKrasList$haigisActivationScore %in%
                                   actScores[actClass=="Strong"] ] = "Strong"


# Codons 12, 13 --------------------------------------------------------



combinedKrasList$class1213 <- rep("Weak",nrow(combinedKrasList))
combinedKrasList$class1213[combinedKrasList$codon %in% c(12,13)] = "Strong"


# Frequency in cbioportal relative to non-hypermutant CRC mut rate --------


# analysis at the amino acid level (so mutcds that give the same amino acid changed are bundled).

krasDataMisNoPoleEdm <- read.csv(file =paste0(reduced_datadir,"cbioKrasDataNoPoleEdm.csv"),
                               header = T)


#get unique patient - mutation pairs
noDupsKrasMisNoPoleEdm <- ddply(krasDataMisNoPoleEdm[,
                                                     c("PatientID",
                                                       "Protein.Change",
                                                       "aaNumber",
                                                       "HGVSc")],
                                .(PatientID,Protein.Change,
                                  aaNumber,HGVSc),nrow)


#annotate with mutCds and mutAA
noDupsKrasMisNoPoleEdm$mutCds <- sapply(noDupsKrasMisNoPoleEdm$HGVSc, function(x){
  mutcds <- paste0("c.",(x %>% strsplit(.,"[.]"))[[1]][3] ) %>% 
    return
})
noDupsKrasMisNoPoleEdm$mutAA <- sapply(noDupsKrasMisNoPoleEdm$Protein.Change, function(x){
  paste0("p.",x) %>% return
})

#only consider mutations that are in the combined Kras list
noDupsKrasMisNoPoleEdmCombList <- subset(noDupsKrasMisNoPoleEdm, 
                                         mutAA %in% combinedKrasList$mutAA)

#construct data frame with mutAA counts
cbioKrasNoPoleMutAA <- data.frame( "Count" = noDupsKrasMisNoPoleEdmCombList$mutAA %>% table)

colnames(cbioKrasNoPoleMutAA ) <- c("mutAA", "count")

cbioKrasNoPoleMutAA$codon <- sapply(as.character(cbioKrasNoPoleMutAA$mutAA), function(x){
  pchange <- (strsplit(x,"[.]"))[[1]][2]
  return(substr(pchange,2,nchar(pchange)-1) %>% as.numeric)
})

cbioKrasNoPoleMutAA <- cbioKrasNoPoleMutAA[order(cbioKrasNoPoleMutAA$codon),]

#get the expected mutation count if mutations were assigned only by 
#mutation rate.
cbioKrasNoPoleMutAA$relMut <- sapply(cbioKrasNoPoleMutAA$mutAA, function(x){
  sum(subset(combinedKrasList,mutAA == x)$relMutRateCrc) %>% return
})


cbioKrasNoPoleMutAA$scaledRelMut <- cbioKrasNoPoleMutAA$relMut/sum(cbioKrasNoPoleMutAA$relMut )

#### Null model for cbiportal data ####
# Take total kras mutations observed and null is that AA muts seen should be
# proportional to mutation rate
cbioKrasNoPoleMutAA$NullCount <- cbioKrasNoPoleMutAA$scaledRelMut*sum(cbioKrasNoPoleMutAA$count)
cbioKrasNoPoleMutAA$countNullRatio <-cbioKrasNoPoleMutAA$count/cbioKrasNoPoleMutAA$NullCount 


assignWeakStrongKmeans <- function(vals){
  kmeanVals <- kmeans(vals,2)
  classVec <- rep("Weak", length(vals))
  if (which.max( kmeanVals$centers)==2){
    #2nd cluster is biggest
    classVec[kmeanVals$cluster==2] = "Strong"
  } else {
    classVec[kmeanVals$cluster==1]  = "Strong"
  }
  return(classVec)
}
cbioKrasNoPoleMutAA$classCbioCnt <- assignWeakStrongKmeans(cbioKrasNoPoleMutAA$countNullRatio)


cbioKrasNoPoleMutAA$index <- 1:nrow(cbioKrasNoPoleMutAA)



#### annotate combined kras list ####
combinedKrasList$classCbioCnt <- sapply(combinedKrasList$mutAA, function(aa){
  
  if (aa %in% cbioKrasNoPoleMutAA$mutAA){
    subset(cbioKrasNoPoleMutAA,mutAA == aa)$classCbioCnt %>% as.character %>% return()
  } else {
    return("Weak") #if AA not seen in cbio data for non pole crcs
  }
  
})






# MSI filtering section: as above with checks for MSI effect ---------------------------------------------------------

#### Filtering: must be known MSS ####


krasDataMisNoPoleEdm <- read.csv(file =paste0(reduced_datadir,"cbioKrasDataNoPoleEdm.csv"),
                                 header = T)

#assessing affect of filtering known MSI status
#note that if delete "" (i.e. samples with unknown status then no change in result)
krasDataMisNoPoleEdmMSS <- subset(krasDataMisNoPoleEdm, !(MSI.Status %in% c("MSI"  ,
                                                                         "MSI-high",
                                                                         "MSI-H",
                                                                         "")))


#get unique patient - mutation pairs
noDupsKrasMisNoPoleEdmMSS <- ddply(krasDataMisNoPoleEdmMSS[,
                                                     c("PatientID",
                                                       "Protein.Change",
                                                       "aaNumber",
                                                       "HGVSc")],
                                .(PatientID,Protein.Change,
                                  aaNumber,HGVSc),nrow)


#annotate with mutCds and mutAA
noDupsKrasMisNoPoleEdmMSS$mutCds <- sapply(noDupsKrasMisNoPoleEdmMSS$HGVSc, function(x){
  mutcds <- paste0("c.",(x %>% strsplit(.,"[.]"))[[1]][3] ) %>% 
    return
})
noDupsKrasMisNoPoleEdmMSS$mutAA <- sapply(noDupsKrasMisNoPoleEdmMSS$Protein.Change, function(x){
  paste0("p.",x) %>% return
})

#only care about mutations that are in the combined Kras list
noDupsKrasMisNoPoleEdmMSSCombList <- subset(noDupsKrasMisNoPoleEdmMSS, 
                                         mutAA %in% combinedKrasList$mutAA)

#construct data frame with mutAA counts
cbioKrasNoPoleMSSMutAA <- data.frame( "Count" = noDupsKrasMisNoPoleEdmMSSCombList$mutAA %>% table)

colnames(cbioKrasNoPoleMSSMutAA ) <- c("mutAA", "count")

cbioKrasNoPoleMSSMutAA$codon <- sapply(as.character(cbioKrasNoPoleMSSMutAA$mutAA), function(x){
  pchange <- (strsplit(x,"[.]"))[[1]][2]
  return(substr(pchange,2,nchar(pchange)-1) %>% as.numeric)
})

cbioKrasNoPoleMSSMutAA <- cbioKrasNoPoleMSSMutAA[order(cbioKrasNoPoleMSSMutAA$codon),]

#get the expected mutation count if mutations were assigned only by 
#mutation rate.
cbioKrasNoPoleMSSMutAA$relMut <- sapply(cbioKrasNoPoleMSSMutAA$mutAA, function(x){
  sum(subset(combinedKrasList,mutAA == x)$relMutRateCrc) %>% return
})


cbioKrasNoPoleMSSMutAA$scaledRelMut <- cbioKrasNoPoleMSSMutAA$relMut/sum(cbioKrasNoPoleMSSMutAA$relMut )


#### Null model for cbiportal data ####
# Take total kras mutations observed and null is that AA muts seen should be
# proportional to mutation rate. cbioKrasNoPoleMSSMutAA
cbioKrasNoPoleMSSMutAA$NullCount <- cbioKrasNoPoleMSSMutAA$scaledRelMut*sum(cbioKrasNoPoleMSSMutAA$count)
cbioKrasNoPoleMSSMutAA$countNullRatio <-cbioKrasNoPoleMSSMutAA$count/cbioKrasNoPoleMSSMutAA$NullCount 



assignWeakStrongKmeans <- function(vals){
  kmeanVals <- kmeans(vals,2)
  classVec <- rep("Weak", length(vals))
  if (which.max( kmeanVals$centers)==2){
    #2nd cluster is biggest
    classVec[kmeanVals$cluster==2] = "Strong"
  } else {
    classVec[kmeanVals$cluster==1]  = "Strong"
  }
  return(classVec)
}
cbioKrasNoPoleMSSMutAA$classCbioCnt <- assignWeakStrongKmeans(cbioKrasNoPoleMSSMutAA$countNullRatio)
strongAAMSScbio <- subset(cbioKrasNoPoleMSSMutAA, classCbioCnt == "Strong")$mutAA %>% as.character()


combinedKrasList$classCbioCntMssStrict = rep("Weak",
                                          nrow(combinedKrasList))
combinedKrasList[combinedKrasList$mutAA %in% strongAAMSScbio,]$classCbioCntMssStrict = "Strong"




# Save outputs ------------------------------------------------------------



if (saveOutput == T){
  
  filenameCombinedKras = paste0(reduced_datadir,"combinedKrasListDriverAnnotated.csv")
  write.csv(combinedKrasList,
            file = filenameCombinedKras,
            row.names =  F)
  
  
}

