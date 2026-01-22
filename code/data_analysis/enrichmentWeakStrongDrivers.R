#Script to assess enrichment of weak/strong drivers in POLE vs non-Pole 
#mutated CRC cancers

#Note we filter cbioportal data by variants in KRAS 
#driver list. Thus, not all mutations observed in 
#ciobio data will be included.

# Preamble ----------------------------------------------------------------
library(readxl)
library(plyr)
library(magrittr)
library(gridExtra)
library(cowplot)

#data downloaded from https://www.cbioportal.org/study/summary?id=coadread_dfci_2016%2Ccoadread_genentech%2Ccoadread_tcga%2Ccoadread_tcga_pub%2Ccoadread_mskcc%2Ccoad_caseccc_2015%2Ccoad_cptac_2019%2Ccoadread_tcga_pan_can_atlas_2018%2Ccoadread_mskresistance_2022%2Ccrc_apc_impact_2020%2Ccrc_dd_2022%2Ccrc_nigerian_2020%2Ccrc_msk_2017%2Crectal_msk_2019%2Crectal_msk_2022%2Cbowel_colitis_msk_2022
#then query kras and pole


reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")
plotDir <- paste0(projectRoot,"images/krasPublicData/")


saveOutput = F

#load in dataframe of kras mutations, annotated with driver info
combinedKrasList <- read.csv(file =paste0(reduced_datadir,"combinedKrasListDriverAnnotated.csv"),
                             header = T)



# Count of mutations in cbio data -----------------------------------------

#### mutation counts in cases without pole edm ####

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

#only care about mutations that are in the combined Kras list
noDupsKrasMisNoPoleEdmCombList <- subset(noDupsKrasMisNoPoleEdm, 
                                         mutCds %in% combinedKrasList$mutCds)


#construct data frame with mutAA counts
cbioKrasNoPoleMutCds <- noDupsKrasMisNoPoleEdmCombList$mutCds %>% 
                                      table

#### mutation counts in cases without pole edm ####


krasDataMisWPoleEdm <- read.csv(file =paste0(reduced_datadir,"cbioKrasDataWPoleEdm.csv"),
                                 header = T)

#get unique patient - mutation pairs
noDupsKrasMisWPoleEdm <- ddply(krasDataMisWPoleEdm[,
                                                     c("PatientID",
                                                       "Protein.Change",
                                                       "aaNumber",
                                                       "HGVSc")],
                                .(PatientID,Protein.Change,
                                  aaNumber,HGVSc),nrow)

#annotate with mutCds and mutAA
noDupsKrasMisWPoleEdm$mutCds <- sapply(noDupsKrasMisWPoleEdm$HGVSc, function(x){
  mutcds <- paste0("c.",(x %>% strsplit(.,"[.]"))[[1]][3] ) %>% 
    return
})

#only care about mutations that are in the combined Kras list
noDupsKrasMisWPoleEdmCombList <- subset(noDupsKrasMisWPoleEdm, 
                                         mutCds %in% combinedKrasList$mutCds)

#construct data frame with mutAA counts
cbioKrasWPoleMutCds <- noDupsKrasMisWPoleEdmCombList$mutCds %>% 
                                      table


# Annotate Kras List with Cbio counts -------------------------------------

combinedKrasList$cbioCntsWPole <- sapply(combinedKrasList$mutCds, function(x){
  if (x %in% names(cbioKrasWPoleMutCds )){
    return(cbioKrasWPoleMutCds[x])
  } else{
    return(0)
  }
})

combinedKrasList$cbioCntsNoPole <- sapply(combinedKrasList$mutCds, function(x){
  if (x %in% names(cbioKrasNoPoleMutCds )){
    return(cbioKrasNoPoleMutCds[x])
  } else{
    return(0)
  }
})


# Fisher test on enrichment -----------------------------------------------


possibleClassVars <- colnames(combinedKrasList)[colnames(combinedKrasList) %>%
                                                  grep("class",.)]

possibleClassVarsExplanation <- c(
  "Bioactivation score",
  "Exon 2 mutations are strong drivers",
  "observed:expected",
  "observed:expected, MSS only"
  
)

#note that while noDupsKrasMisWPoleEdmCombList has 31 entries
#not all of these mutations are considered in combinedKRAS list of drivers
#hence have less than 31 in sum of first row in contingency table
contingencyTableClassVars <- lapply(possibleClassVars, function(x){
  numWPoleWeak <- combinedKrasList$cbioCntsWPole[combinedKrasList[,x]=="Weak"] %>% sum
  numWPoleStrong <- combinedKrasList$cbioCntsWPole[combinedKrasList[,x]=="Strong"] %>% sum
  
  numNoPoleWeak <- combinedKrasList$cbioCntsNoPole[combinedKrasList[,x]=="Weak"] %>% sum
  numNoPoleStrong <- combinedKrasList$cbioCntsNoPole[combinedKrasList[,x]=="Strong"] %>% sum
  
  contTable <- rbind(c(numWPoleWeak,numWPoleStrong),
                     c(numNoPoleWeak,numNoPoleStrong))
  
  rownames(contTable) <- c("Pole","Non-Pole")
  colnames(contTable) <- c("Weak driver", "Strong driver")

  return(contTable)
})

fisherTest <- lapply(1:length(contingencyTableClassVars), function(k){
  fisher.test(contingencyTableClassVars[[k]]) %>% return()
})
#pvalues all <.05



contingencyTableGridPlot <- lapply(1:length(contingencyTableClassVars), function(k){
  
  classTitle <- ggdraw() + 
    draw_label(paste0(possibleClassVarsExplanation[k],
                      "\nOR =",as.character(fisherTest[[k]]$estimate),
                      "\np-value =",as.character(fisherTest[[k]]$p.value)),
               fontface='bold')
  
  grid.arrange(classTitle,
               tableGrob(contingencyTableClassVars[[k]]),ncol = 2) %>% return()
  # textGrob("Daily QC: Blue",gp=gpar(fontsize=20,font=3))
  
})

contingencyPlotAll <- grid.arrange(contingencyTableGridPlot[[1]],
                     contingencyTableGridPlot[[2]],
                     contingencyTableGridPlot[[3]],
                     contingencyTableGridPlot[[4]],
                     nrow =4)





# Location of tumour ------------------------------------------------------
#possible columns in cbio with site name. Just use Primary.Tumor.Site.
siteCNames <- c("Primary.Tumor.Site"    ,
                "Tumor Site"    ,
                "International Classification of Diseases for Oncology, Third Edition ICD-O-3 Site Code",
                "Patient Primary.Tumor.Site",
                "Tumor Disease Anatomic Site" ,
                "Tissue Site" ,
                "Tissue site of derivation")

distalTerms <- c("3 - left colon" , "Rectum", "Sigmoid Colon","Rectum","Descending Colon","Left","4 - rectum", "Left Colon",
                 "Splenic Flexure","Rectosigmoid" ,"Descending","Rectosigmoid Colon" ,"Sigmoid" )
proximalTerms <- c("1 - right colon", "Ascending Colon","Cecum","Transverse Colon","Right","2 - transverse colon",
                   "R Colon","Hepatix Flexure" ,"Ascending colon"  )


## annotate POLE mutated samples
noDupsKrasMisWPoleEdmSite <- ddply(krasDataMisWPoleEdm[,c("PatientID","Protein.Change","aaNumber","Primary.Tumor.Site")],
                                   .(PatientID,`Protein.Change`,aaNumber,`Primary.Tumor.Site`),nrow)
noDupsKrasMisWPoleEdmWSite <-  subset(noDupsKrasMisWPoleEdmSite,`Primary.Tumor.Site`%in% c(distalTerms,proximalTerms))

noDupsKrasMisWPoleEdmWSite$distProx <- rep(NA,nrow(noDupsKrasMisWPoleEdmWSite))
noDupsKrasMisWPoleEdmWSite$distProx[noDupsKrasMisWPoleEdmWSite$`Primary.Tumor.Site` %in% distalTerms] = "distal"
noDupsKrasMisWPoleEdmWSite$distProx[noDupsKrasMisWPoleEdmWSite$`Primary.Tumor.Site` %in% proximalTerms] = "proximal"

nMutsWPoleProximal <- sum(noDupsKrasMisWPoleEdmWSite$distProx == "proximal")
nMutsWPoleDistal <- sum(noDupsKrasMisWPoleEdmWSite$distProx == "distal")


aaNumsWPoleDistal <- subset(noDupsKrasMisWPoleEdmWSite,distProx == "distal")$aaNumber
aaNumsWPoleProximal <- subset(noDupsKrasMisWPoleEdmWSite,distProx == "proximal")$aaNumber

percMutsWPoleDistal <- (aaNumsWPoleDistal %>% table)/nMutsWPoleDistal
percMutsWPoleProximal <- (aaNumsWPoleProximal %>% table)/nMutsWPoleProximal

## annotate non-POLE mutated samples
noDupsKrasMisNoPoleEdmSite <- ddply(krasDataMisNoPoleEdm[,c("PatientID","Protein.Change","aaNumber","Primary.Tumor.Site")],
                                    .(PatientID,`Protein.Change`,aaNumber,`Primary.Tumor.Site`),nrow)
noDupsKrasMisNoPoleEdmWSite <-  subset(noDupsKrasMisNoPoleEdmSite,`Primary.Tumor.Site` %in% c(distalTerms,proximalTerms))

noDupsKrasMisNoPoleEdmWSite$distProx <- rep(NA,nrow(noDupsKrasMisNoPoleEdmWSite))
noDupsKrasMisNoPoleEdmWSite$distProx[noDupsKrasMisNoPoleEdmWSite$`Primary.Tumor.Site` %in% distalTerms] = "distal"
noDupsKrasMisNoPoleEdmWSite$distProx[noDupsKrasMisNoPoleEdmWSite$`Primary.Tumor.Site` %in% proximalTerms] = "proximal"

nMutsNoPoleProximal <- sum(noDupsKrasMisNoPoleEdmWSite$distProx == "proximal")
nMutsNoPoleDistal <- sum(noDupsKrasMisNoPoleEdmWSite$distProx == "distal")

aaNumsNoPoleDistal <- subset(noDupsKrasMisNoPoleEdmWSite,distProx == "distal")$aaNumber
aaNumsNoPoleProximal <- subset(noDupsKrasMisNoPoleEdmWSite,distProx == "proximal")$aaNumber

percMutsNoPoleDistal <- (aaNumsNoPoleDistal %>% table)/nMutsNoPoleDistal
percMutsNoPoleProximal <- (aaNumsNoPoleProximal  %>% table)/nMutsNoPoleProximal


#### test enrichment of codon 12/13 muts by distal proximal ####
contTableExon2Distal <- rbind( c(sum ( !(aaNumsWPoleDistal  %in% c(12,13))) , 
                                   sum (aaNumsWPoleDistal  %in% c(12,13))),
                               c( sum ( !(aaNumsNoPoleDistal  %in% c(12,13))), 
                                   sum (aaNumsNoPoleDistal  %in% c(12,13))) )

rownames(contTableExon2Distal) <- c("Pole","Non-Pole")
colnames(contTableExon2Distal) <- c("Weak driver", "Strong driver")

fishTest1213Distal <- fisher.test(contTableExon2Distal)

perc1213WPoleDistal <- sum (aaNumsWPoleDistal  %in% c(12,13))/nMutsWPoleDistal
perc1213NoPoleDistal <- sum (aaNumsNoPoleDistal  %in% c(12,13))/nMutsNoPoleDistal

contTableExon2Proximal <- rbind( c(sum ( !(aaNumsWPoleProximal %in% c(12,13))),
                                  sum (aaNumsWPoleProximal %in% c(12,13))),
                               c(  sum ( !(aaNumsNoPoleProximal  %in% c(12,13)))  , 
                                 sum (aaNumsNoPoleProximal  %in% c(12,13))))

rownames(contTableExon2Proximal) <- c("Pole","Non-Pole")
colnames(contTableExon2Proximal) <- c("Weak driver", "Strong driver")
fishTest1213Proximal <- fisher.test(contTableExon2Proximal)
perc1213WPoleProximal <- sum (aaNumsWPoleProximal   %in% c(12,13))/nMutsWPoleProximal 
perc1213NoPoleProximal  <- sum (aaNumsNoPoleProximal   %in% c(12,13))/nMutsNoPoleProximal 

locExplan<- c(
  "Distal CRC\n",
  "Proximal CRC\n"
)

fishersTestsLoc = list(fishTest1213Distal,
                       fishTest1213Proximal)

contTableLocs <- list(contTableExon2Distal,contTableExon2Proximal)
contingencyTableGridLocPlot <- lapply(1:2, function(k){
  

  
  classTitle <- ggdraw() + 
    draw_label(paste0(locExplan[k],
                             "Exon 2 mutations are strong drivers",
                      "\nOR =",as.character(fishersTestsLoc[[k]]$estimate),
                      "\np-value =",as.character(fishersTestsLoc[[k]]$p.value)),
               fontface='bold')
  
  grid.arrange(classTitle,
               tableGrob(contTableLocs[[k]]) ,ncol = 2) %>% return()
  # textGrob("Daily QC: Blue",gp=gpar(fontsize=20,font=3))
  
})
grid.arrange(contingencyTableGridLocPlot[[1]],
             contingencyTableGridLocPlot[[2]],
             nrow =2)


contingencyPlotAll <- grid.arrange(contingencyTableGridPlot[[1]],
                                   contingencyTableGridPlot[[2]],
                                   contingencyTableGridPlot[[3]],
                                   contingencyTableGridPlot[[4]],
                                   contingencyTableGridLocPlot[[1]],
                                   contingencyTableGridLocPlot[[2]],
                                   nrow =6)



#### save contingency table plots ####

if (saveOutput == T){
  pltfilename <- paste0(plotDir,"contingencyPlotAll.pdf")
  ggsave2(filename = pltfilename,contingencyPlotAll,
          height= 20,
          width = 20,
          units = "cm")
  
  filename <- paste0(reduced_datadir,"contigencyTables.rds")
  saveRDS(contingencyTableClassVars, file = filename)
  
  
  
  
  
}
