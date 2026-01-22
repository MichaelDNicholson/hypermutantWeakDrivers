# annotate KRAS cbioportal data for whether samples
# were POLE WT/mutant

#Bottom of script output summary of studies used

library(magrittr)
library(data.table)
library(plyr)


# Preamble ----------------------------------------------------------------

#data downloaded from https://www.cbioportal.org/study/summary?id=coadread_dfci_2016%2Ccoadread_genentech%2Ccoadread_tcga%2Ccoadread_tcga_pub%2Ccoadread_mskcc%2Ccoad_caseccc_2015%2Ccoad_cptac_2019%2Ccoadread_tcga_pan_can_atlas_2018%2Ccoadread_mskresistance_2022%2Ccrc_apc_impact_2020%2Ccrc_dd_2022%2Ccrc_nigerian_2020%2Ccrc_msk_2017%2Crectal_msk_2019%2Crectal_msk_2022%2Cbowel_colitis_msk_2022
#then query kras and pole
cbioportalDataDir <- paste0(projectRoot,"data/cbioportal/")
reduced_datadir <- paste0(projectRoot,"reduced_data/")
plotdir = paste0(projectRoot,"images/krasPublicData/")
fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveplots = F
saveOutput = F



# Functions ---------------------------------------------------------------

getAAnumberPchange <- function(pchange){
  if (pchange==""){
    return(NA)
  } else {
    return(substr(pchange,2,nchar(pchange)-1) %>% as.numeric())
  }
}

# Annotate samples with pole exonuclease domain mutation ----------------------------------------------------------------
# use clinical data to link samples with pole & kras muts

clinicalData <- as.data.frame(fread(paste0(cbioportalDataDir,
                                           "cbioportal_02052023_CRC_combined_study_clinical_data.tsv")))
krasData <- as.data.frame(fread(paste0(cbioportalDataDir,
                                       "cbioportal_02052023_CRC_KRAS_mutdata.tsv")))
poleData <- as.data.frame(fread(paste0(cbioportalDataDir,
                                       "cbioportal_02052023_CRC_POLE_mutdata.tsv")))

#list of pathogenic Pole mutations
pathoPoleMuts <- read.csv(paste0(reduced_datadir,"pathogenicPoleList.csv"))[,2]


poleDataMis <- subset(poleData, `Mutation Type` == "Missense_Mutation" &
                        !grepl("del",`Protein Change`) &
                        !grepl("del",HGVSc))

# exonuclease domain aa number https://academic.oup.com/hmg/article/22/14/2820/752894
poleDataMis$Edm <- poleDataMis$`Protein Change` %in% pathoPoleMuts
samplesWPoleEdm <- subset(poleDataMis,Edm == T &
                            `Protein Change` %in%  pathoPoleMuts)$`Sample ID`


# Annotate KRAS mutations by AA change & Pole status ----------------------------------------------------------------

krasDataMis <- subset(krasData, `Mutation Type` == "Missense_Mutation" &
                        !grepl("del",`Protein Change`) &
                        !grepl("del",HGVSc)) #extra filter for indel

#Interested in unique patient-mutation combinations
# annotate KRAS muts with Patient ID
krasDataMis$PatientID <- sapply(krasDataMis$`Sample ID`, function(x){
  subset(clinicalData, `Sample ID` ==x)$`Patient ID` %>% 
    unique %>% 
    return
})

krasDataMis$aaNumber <- sapply(krasDataMis$`Protein Change`, getAAnumberPchange)


krasDataMisWPoleEdm <- subset(krasDataMis, `Sample ID` %in% samplesWPoleEdm)
krasDataMisNoPoleEdm <- subset(krasDataMis, !(`Sample ID` %in% samplesWPoleEdm))


if (saveOutput == T){
  
  
  #kras mutations with Pole edm
  #needed for enrichment of weak/strong in 
  filenameKrasWPole = paste0(reduced_datadir,"cbioKrasDataWPoleEdm.csv")
  write.csv(krasDataMisWPoleEdm ,
            file = filenameKrasWPole ,
            row.names =  F)
  
  #kras mutations seen in samples without POLE edm
  #can be used for classifying weak/strong driver (ob:ep analysis)
  
  filenameKrasNoPole = paste0(reduced_datadir,"cbioKrasDataNoPoleEdm.csv")
  write.csv(krasDataMisNoPoleEdm,
            file = filenameKrasNoPole,
            row.names =  F)
}




# Studies used ------------------------------------------------------------
cBioStudySampUni <- clinicalData[,c("Study ID", "Sample ID")] %>% unique()
cBioStudySampUnTable <- cBioStudySampUni$`Study ID` %>% table

names(cBioStudySampUnTable) <- sapply(names(cBioStudySampUnTable), function(n){
  paste("cBioPortal: ",n)
})

krasData100kgp <- read.csv(file =paste0(projectRoot,
                                        "data/100KGP/krasCodonsMutStatusCounts_v8_08012026.csv"),
                           header = T)[,-1]

n100KgpSampsUsed <- krasData100kgp[,c("MSS","MSI","MSI.CIN")] %>% rowSums(na.rm=T) %>% sum()

combSummaryStudySamples <- data.frame(cBioStudySampUnTable)
colnames(combSummaryStudySamples) <- c("Study","Samples")
combSummaryStudySamples <- rbind(combSummaryStudySamples, data.frame("Study" = c("100,000 genomes project V8 (only with oncogenic KRAS)"),
                                                                       "Samples" = n100KgpSampsUsed  ))



if (saveOutput == T){
  
  

  
  filenameStudySummary = paste0(reduced_datadir,"studySummaryTable.csv")
  write.csv(combSummaryStudySamples,
            file = filenameStudySummary,
            row.names =  F)
}
