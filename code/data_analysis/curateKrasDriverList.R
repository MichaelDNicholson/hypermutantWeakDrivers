# Script to curate list of kras driver variants
# Combine driver list from intogen and clinical genetics input
# Require variables mutCds mutAA triContext mutType codon boostDM
# Result: combinedKrasList, variants either intogen CRC driver, oncogenic or 
# likely oncogenic in OncoKB (Jan 2025) + D57N (observed in 2 cbio pole samples
# +1 non-pole)

# Then annotate with relative mutation rate for each mutation
# under non-hypermutant spectra and pole-mutant

# Preamble ----------------------------------------------------------------

library(magrittr)
library("BSgenome.Hsapiens.UCSC.hg38")
library(readxl)



reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")
miscelDir <- paste(projectRoot,"data/miscellaneous/",sep = "")
saveOutput <- F


pyrBases <- c("C","T")
watsonCrickKey = data.frame("A" = "T", "C" = "G", "T" = "A", "G" = "C")

#use cosmic mutation list to annotate entries.
cosmicKrasList <-  read.csv(paste0(miscelDir,"V94_38_MUTANTCENSUS_KRAS.csv")) 


#### Load in driver lists ###
intogenKras <- read.table(paste0(miscelDir,"IntOGen-Distribution-KRAS-COREAD.tsv"),header = T)
intogenKrasDriver <- subset(intogenKras,Driver == "Driver") #filter by driver classification according to boostDM
customKrasList = read_excel(paste0(miscelDir,"custom_KRAS_driverlist.xlsx"))

#### Load in representative mutation spectra and trinucleotide counts ###
repCrcPoleSpectra <- read.csv(paste0(reduced_datadir,
                                     "pcawgRepPoleSpectra.csv"),
                              header = T)
repCrcSpectra <- read.csv(paste0(reduced_datadir,
                                 "pcawgRepNonHypermutSpectra.csv"),
                          header = T)
repMsiSpectra <- read.csv(paste0(reduced_datadir,
                                 "repMsiSpectra.csv"),
                          header = T)

triCountsGenome <- read.csv(paste0(reduced_datadir,
                                   "tricntsHg38PriAssembly.csv"),
                            header = T)





# Functions ---------------------------------------------------------------

get_mutType = function(mutpos,mutCds){
  #this relies on knowledge that coding strand for kras is on reverse strand
  # so if pyrimidine is in mutCds
  # tricontext is rev(compliment(forward[pos-1:pos+1]))
  # if pyr not in mutcds
  # tricontext is forward[pos-1:pos+1]
  # returns of form X[M>N]Y
  
  #example input
  # mutpos <- intogenKrasDriver$mutPos[5]
  # mutCds <- intogenKrasDriver$mutCds[5]
  
  ineqLoc <- gregexpr(">",mutCds)[[1]][1]
  baseChange <- substr(mutCds,ineqLoc-1,ineqLoc+1)
  refBase <- substr(baseChange,1,1)
  if (refBase %in% pyrBases ){
    # pyr is on reverse (kras coding on reverse)
    #depends on refBase is forward or reverse
    
    mutType <- paste0(watsonCrickKey[ Hsapiens$chr12[mutpos+1] %>% as.character()],
                      "[",
                      baseChange,
                      "]",
                      watsonCrickKey[ Hsapiens$chr12[mutpos-1] %>% as.character()])
  } else {
    #pyr is on forward strand
    baseChangeSplit = strsplit(baseChange,">") %>% unlist
    baseChangePyr <- paste0(watsonCrickKey[baseChangeSplit[1]],
                            ">",
                            watsonCrickKey[baseChangeSplit[2]])
    
    mutType <- paste0(Hsapiens$chr12[mutpos-1],
                      "[",
                      baseChangePyr,
                      "]",
                      Hsapiens$chr12[mutpos+1])
    
  }
  return(mutType)
}

# Intogen list annotate ---------------------------------------------------

# Use cosmic mutation data for annotation

#for intogen need to link position to cds. Then use cosmic for mut_aa
intogenKrasDriver$mutCds <- sapply(intogenKrasDriver$Mutation..GRCh38., function(mut){
  mutpos <- strsplit(mut,"[:]")[[1]][[2]] %>% as.numeric
  baseChange <- strsplit(mut,"[:]")[[1]][[3]]
  baseChangeSplit <- baseChange %>% strsplit(.,">") %>% unlist
  compBaseChange = paste0(watsonCrickKey[baseChangeSplit[1]],
                          ">",
                          watsonCrickKey[baseChangeSplit[2]])
  #for kras: coding strand is on reverse strand, hence
  #reference base for mut_cds taken from reverse (i.e. why we use compBaseChange)
  subset(cosmicKrasList,grepl(mutpos,MUTATION_GENOME_POSITION) &
           grepl(compBaseChange,MUTATION_CDS)  )$MUTATION_CDS %>% 
    unique %>% 
    return()
})

intogenKrasDriver$mutAA <- sapply(intogenKrasDriver$mutCds, function(x) {
  subset(cosmicKrasList, MUTATION_CDS == x)$MUTATION_AA %>% unique %>% return
})

intogenKrasDriver$mutPos <- sapply(intogenKrasDriver$Mutation..GRCh38., function(mut){
  mutpos <- strsplit(mut,"[:]")[[1]][[2]] %>% as.numeric
})

intogenKrasDriver$mutType <- sapply(1:nrow(intogenKrasDriver), function(k){
  return(get_mutType(intogenKrasDriver$mutPos[k],
                     intogenKrasDriver$mutCds[k]))
})

intogenKrasDriver$triContext <- sapply(intogenKrasDriver$mutType , function(x){
  paste0(substr(x,1,1),
         substr(x,3,3),
         substr(x,7,7)) %>% 
    return
})

#is mutation classified as driver in the intogen list
intogenKrasDriver$boostDM <- rep("TRUE",nrow(intogenKrasDriver))

# custom list annotate ------------------------------------------------

customKrasList$mutCds = sapply(customKrasList$CDS,function(x) 
  paste("c.",paste( strsplit(x," ") %>% unlist,collapse =  ""),sep ="") %>% return)

customKrasList$mutAA = sapply(customKrasList$aa_change,function(x) paste0("p.",x))

customKrasList$mutPos <- sapply(customKrasList$mutCds , function(x){
  relMutGenPos <-  subset(cosmicKrasList,MUTATION_CDS == x)$MUTATION_GENOME_POSITION %>% unique()
  ( (relMutGenPos %>% strsplit(.,":"))[[1]][2] %>% strsplit(.,"-") )[[1]][1] %>% 
    as.numeric %>% 
    return()
})


customKrasList$mutType <- sapply(1:nrow(customKrasList), function(k){
  return(get_mutType(customKrasList$mutPos[k],
                     customKrasList$mutCds[k]))
})

customKrasList$triContext <- sapply(customKrasList$mutType , function(x){
  paste0(substr(x,1,1),
         substr(x,3,3),
         substr(x,7,7)) %>% 
    return
})

#is mutation classified as driver in the intogen list
customKrasList$boostDM <- customKrasList$mutCds %in% intogenKrasDriver$mutCds

# Combine Kras driver lists and save --------------------------------------
commonCols <- c("mutCds",
                "mutAA",
                "triContext",
                "mutType",
                "boostDM")

combinedKrasList <- rbind(intogenKrasDriver[,commonCols],
                          customKrasList[,commonCols]) %>% 
  unique

combinedKrasList$codon <- sapply(combinedKrasList$mutAA, function(x){
  substr(x,4,nchar(x)-1) %>% as.numeric %>% return
})

combinedKrasList <- combinedKrasList[order(combinedKrasList$codon),]



# Relative mutation rate --------------------------------------------------
# When mutation occurs it is of given mutation type X[A>B]Y with probability
# proportional to spectra values (other factors common to both mutation subtypse, e.g. 
# number of matching tricontext.



#### Get folded tricounts genome wide ###
#record only tricontext with pyr central base
#record 5' and 3' bases

triCountsGenome$foldedContext <- sapply(triCountsGenome$Context, function(x){
  centralBase <- substr(x,2,2)
  if (centralBase %in% pyrBases ){
    return(x)
  } else { 
    complementContext <-  paste0(watsonCrickKey[substr(x,3,3)],
                                 watsonCrickKey[substr(x,2,2)],
                                 watsonCrickKey[substr(x,1,1)])
    return(complementContext)
  }
})

triCountsFoldedGenome <- sapply(triCountsGenome$foldedContext %>% unique, function(x){
  sum(subset(triCountsGenome, foldedContext == x)$Count) %>% return
})

#### Add relative mutation rate for different spectra ####

#relative mutation rate assuming POLE mutational process
combinedKrasList$relMutRateCrcPole = sapply(1:nrow(combinedKrasList), function(k){
  tri <-  combinedKrasList[k,"triContext"]
  mutclass <-  combinedKrasList[k,"mutType"]
  #prob mut falls on in mut class (tri and base change)
  probMutType <-   subset(repCrcPoleSpectra, MutationType == mutclass)$SpectraVals 
#   return( probMutType )
# })
  #Given of mut type, probability it occurs at that location.
  probAtLoc <- (1/triCountsFoldedGenome[tri]) 
  return( probMutType *probAtLoc )
})

#relative mutation rate assuming non-hypermutant CRC mutational processes
combinedKrasList$relMutRateCrc = sapply(1:nrow(combinedKrasList), function(k){
  tri <-  combinedKrasList[k,"triContext"]
  mutclass <-  combinedKrasList[k,"mutType"]
  #prob mut falls on in mut class (tri and base change)
  probMutType <-   subset(repCrcSpectra, MutationType == mutclass)$SpectraVals 
#   return( probMutType )
# })
    probAtLoc <-  (1/triCountsFoldedGenome[tri]) 
  return( probMutType *probAtLoc )
})


# For MSI --------------------------------------------------

combinedKrasList$relMutRateMSI= sapply(1:nrow(combinedKrasList), function(k){
  tri <-  combinedKrasList[k,"triContext"]
  mutclass <-  combinedKrasList[k,"mutType"]
  #prob mut falls on in mut class (tri and base change)
  probMutType <-   subset(repMsiSigdf , MutationType == mutclass)$SpectraVals 
  #   return( probMutType )
  # })
  probAtLoc <-  (1/triCountsFoldedGenome[tri]) 
  return( probMutType *probAtLoc )
})


if (saveOutput == T){
  
  filenameCombinedKras = paste0(reduced_datadir,"combinedKrasList.csv")
  write.csv(combinedKrasList,
            file = filenameCombinedKras,
            row.names =  F)
}


