
# Script to annotate Kras mutations as weak/strong drivers
# Will add columns with weak/strong under varied criteria.

# Preamble ----------------------------------------------------------------

library(magrittr)



reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")

dataDir <- paste0(projectRoot,"data/")
#load in dataframe of kras mutations, annotated with driver info
combinedKrasList <- read.csv(file =paste0(reduced_datadir,"combinedKrasListDriverAnnotated.csv"),
                             header = T)

fnCleanAnnotateKrasData <- function(krasDataSub){
  noDupsKras  <- ddply(krasDataSub[,
                                                             c("PatientID",
                                                               "Protein.Change",
                                                               "aaNumber",
                                                               "HGVSc")],
                                     .(PatientID,Protein.Change,
                                       aaNumber,HGVSc),nrow)
  
  
  
  noDupsKras$mutAA <- sapply( noDupsKras$Protein.Change, function(x){
    paste0("p.",x) %>% return
  })
  
  noDupsKras$codon <- sapply(as.character( noDupsKras$mutAA), function(x){
    pchange <- (strsplit(x,"[.]"))[[1]][2]
    return(substr(pchange,2,nchar(pchange)-1) %>% as.numeric)
  })
  
  #annotate with mutCds and mutAA
  noDupsKras$mutCds <- sapply(noDupsKras$HGVSc, function(x){
    mutcds <- paste0("c.",(x %>% strsplit(.,"[.]"))[[1]][3] ) %>% 
      return
  })
  
  #only care about mutations that are in the combined Kras list
  noDupsKrasCombList <- subset(noDupsKras, 
                                           mutCds %in% combinedKrasList$mutCds)
  
  
  return( noDupsKrasCombList)
}

#### Filtering: MSI and no POLE mut ####


krasDataMisNoPoleEdm <- read.csv(file =paste0(reduced_datadir,"cbioKrasDataNoPoleEdm.csv"),
                                 header = T)
#load in dataframe of kras mutations, annotated with driver info
combinedKrasList <- read.csv(file =paste0(reduced_datadir,"combinedKrasListDriverAnnotated.csv"),
                             header = T)

#assessing affect of filtering known MSI status
#note that if delete "" (i.e. samples with unknown status then no change in result)
krasDataMisNoPoleEdmMSI <- subset(krasDataMisNoPoleEdm, (MSI.Status %in% c("MSI"  ,
                                                                            "MSI-high",
                                                                            "MSI-H")) |
                                    MSI.Type  == "Instable")

krasDataMisNoPoleEdmNoMSI <- subset(krasDataMisNoPoleEdm, !(MSI.Status %in% c("MSI"  ,
                                                                           "MSI-high",
                                                                           "MSI-H")) &
                                    MSI.Type  != "Instable")

noDupsKrasMisNoPoleEdmMSI <- fnCleanAnnotateKrasData(krasDataMisNoPoleEdmMSI )
noDupsKrasMisNoPoleEdmNoMSI <- fnCleanAnnotateKrasData(krasDataMisNoPoleEdmNoMSI  )


msiCodonTable <- noDupsKrasMisNoPoleEdmMSI$codon %>% table
noMsiCodonTable <- noDupsKrasMisNoPoleEdmNoMSI$codon %>% table

msiCodon1213NumCbio <- msiCodonTable[c("12", "13")] %>% sum
msiCodonNot1213NumCbio <- msiCodonTable[!(names(msiCodonTable) %in% c("12", "13"))] %>% 
  sum

nonhypCodon1213NumCbio <- noMsiCodonTable[c("12", "13")] %>% sum
nonhypCodonNot1213NumCbio <- noMsiCodonTable[!(names(noMsiCodonTable) %in% c("12", "13"))] %>% 
  sum




# Add in 100KGP data ------------------------------------------------------
krasData100kgp <- read.csv(file =paste0(dataDir,
                                        "100KGP/krasCodonsMutStatusCounts_v8_08012026.csv"),
                                 header = T)[,-1]
msiCodon1213Num100kgp <-krasData100kgp[krasData100kgp$codon %in% c(12,13),c("MSI","MSI.CIN")] %>% 
  rowSums() %>% sum() 
msiCodonNot1213Num100kgp <-krasData100kgp[!(krasData100kgp$codon %in% c(12,13)),
                                          c("MSI","MSI.CIN")] %>% 
  rowSums(., na.rm= T) %>% sum() 

mssCodon1213Num100kgp<-krasData100kgp[krasData100kgp$codon %in% c(12,13),"MSS"] %>% 
 sum() 
mssCodonNot1213Num100kgp <-krasData100kgp[!(krasData100kgp$codon %in% c(12,13)),"MSS"] %>% 
 sum() 
percMsiNotCodon1213 <- ((msiCodonNot1213NumCbio+msiCodonNot1213Num100kgp)/
  (msiCodonNot1213NumCbio+msiCodonNot1213Num100kgp + 
     msiCodon1213NumCbio+msiCodon1213Num100kgp) ) %>% print
            
fisher.test(rbind(c(msiCodon1213NumCbio+msiCodon1213Num100kgp,
                    msiCodonNot1213NumCbio+msiCodonNot1213Num100kgp),
                  c(nonhypCodon1213NumCbio+mssCodon1213Num100kgp,
                    nonhypCodonNot1213NumCbio+mssCodonNot1213Num100kgp))) %>% print
