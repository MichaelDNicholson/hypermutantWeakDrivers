#Script to construct representative CRC mutational spectra based on
#data from PCAWG
#Aims & structure of script:
# - mixed signature for POLE mutant CRCs
# - signature for non-hypermutant POLE CRCs
# - comparison of signatures

# Preamble ----------------------------------------------------------------

library(magrittr)
library(MutationalPatterns)
library(cowplot)
library(ggplot2)




pcawgDataDir <-  paste0(projectRoot ,
                       "data/pcawg_data/")
cosmicDataDir <-   paste0(projectRoot ,
                          "data/cosmic/")
reduced_datadir = paste(projectRoot,"reduced_data/",sep = "")
saveplots = F
save_reduceddata = F
plotdir =  paste0(projectRoot,"images/pcawg_trinuc_plots/")
todaysDate <- gsub("-","",Sys.Date()) #for filenames

# POLE CRC representative signature ----------------------------------------------------------------

#### Find pcawg samples with POLE mutation ####
#Mutation data for CRC pcawg samples with a POLE mutation. Downloaded from:
# https://dcc.icgc.org/search?filters=%7B%22donor%22:%7B%22studies%22:%7B%22is%22:%5B%22PCAWG%22%5D%7D,%22primarySite%22:%7B%22is%22:%5B%22Colorectal%22%5D%7D%7D,%22gene%22:%7B%22id%22:%7B%22is%22:%5B%22ENSG00000177084%22%5D%7D%7D%7D&donors=%7B%22from%22:1%7D
mutdata = read.table(file = paste0(pcawgDataDir,
                                   'icgc_pcawg_crc_pole_sept21_simple_somatic_mutation.open.tsv.gz'),
                     sep = '\t', header = TRUE)
pole_ens = "ENSG00000177084"
mutdata_pole = subset(mutdata,gene_affected == pole_ens)

#pathogenic POLE mutations
pathoPoleMuts=read.csv(paste0(reduced_datadir,"pathogenicPoleList.csv"),
           row.names = NULL)[,2]

mutdata_pole_patho = subset(mutdata_pole, aa_mutation %in% pathoPoleMuts)
pole_specids = (mutdata_pole_patho$icgc_specimen_id %>% unique())

#### Get signature exposures in POLE PCAWG samples ####

#downloaded from https://www.synapse.org/#!Synapse:syn11804040
sigExposures = read.csv(paste0(pcawgDataDir,
                                  "PCAWG_sigProfiler_SBS_signatures_in_samples.csv"),
                             header = T)

#sigprof samples that correspond to samples with POLE
sigExposuresPoleSamples = subset(sigExposures  ,Sample.Names %in% 
                            pole_specids )

#10a, 10b, 28 are pole associated in cosmic v3
sbsnamesPole <- c("SBS10a", "SBS10b", "SBS28")
sigExposuresPoleSamplesPoleProcess = sigExposuresPoleSamples[,sbsnamesPole]
#normalise signature exposures
sigExposureProbsPole = lapply(1:nrow(sigExposuresPoleSamplesPoleProcess), function(k){
  return(sigExposuresPoleSamplesPoleProcess[k,]/sum(sigExposuresPoleSamplesPoleProcess[k,]))
}) %>% do.call(rbind,.)
#when mut occurs probability it occurs due to each process (10a, 10b, 28)
sigExposureProbsAvPole = colMeans(sigExposureProbsPole)


#### Construct representative POLE signature ####

cosmicSbsVals = read.table(paste0(cosmicDataDir,"COSMIC_v3_SBS_GRCh37.txt"),
                             header = T)

#cosmic v3 contains sbs 84/85 which aren't in the pcawg calls
#but cosmic v2 doesn't contain many pcawg columns- so use v3 but don't use 84/85
cosmicSbsValsNo8485 = subset(cosmicSbsVals,select = -c(SBS84,SBS85))  

repPoleSpectra = sapply(cosmicSbsValsNo8485$Type, function(type){
  
  #get prob mut is of type given mut from pole process x
  perclass_reconstruct = sapply(sbsnamesPole, function(x){
    sigExposureProbsAvPole[x]*cosmicSbsValsNo8485[cosmicSbsValsNo8485$Type==type,x] %>% return()
  })
  return(sum(perclass_reconstruct %>% as.numeric))
})
repPoleSpectraDf = data.frame("MutationType" = names(repPoleSpectra),
                              "SpectraVals" = as.numeric(repPoleSpectra) )

repPoleSpectraVec = as.numeric(repPoleSpectra) %>% as.matrix(.,col = 1)
rownames(repPoleSpectraVec) = repPoleSpectraDf$MutationType
colnames(repPoleSpectraVec) = ("POLE Mutant CRC Spectra")


#Fig 5B lower
plotRepPoleSpectra = plot_96_profile(repPoleSpectraVec,ymax = 0.4)+
  theme(strip.text.y = element_blank(),
        text= element_text(family = "arial",
                           size= 6)
  )


if (saveplots == T){
  pltfilename = paste(plotdir,"repPOLESpectra_",
                todaysDate,".pdf",sep = "")
  
  #hard to save to sizes want in paper
  # save with aspect ratio want and shrink down in inkscape

  ggplot2::ggsave(filename =   pltfilename , 
                  plot =  plotRepPoleSpectra , 
                  device = cairo_pdf, 
                  dpi = 200, 
                  width = 7,
                  height = 2, 
                  units = "in")

  
  }

if (save_reduceddata == T){
  filenameRepPoleSpectraCsv = paste(reduced_datadir,"pcawgRepPoleSpectra.csv",sep="")
  write.csv(repPoleSpectraDf,
            file =filenameRepPoleSpectraCsv,
            row.names =  F)
}



# Non-hypermutant CRC representative signature ----------------------------------------------------------------

#available here: https://www.synapse.org/#!Synapse:syn11726620
pcawgTriCounts =read.csv(unz(paste0(pcawgDataDir,
                                    "WGS_PCAWG_2018_02_09.zip"),
                                    "WGS_PCAWG.96.csv"),header = T)
pcawgCrcTriCounts <- cbind(pcawgTriCounts[,c(1,2)],
                           pcawgTriCounts[,grepl("ColoRect",colnames(pcawgTriCounts))])

##### Exclude hypermutant tumours ####
#define hypermutant tumours based on k-means clustering of total count

totalTriCounts <-  colSums(pcawgCrcTriCounts[,3:ncol(pcawgCrcTriCounts)])
TriCountsSpec <- names(totalTriCounts) %>% gsub("ColoRect.AdenoCA..","",.)

#exclude already identified pole edm containing crcs
totalTriCountsNoPole <- totalTriCounts[!(TriCountsSpec %in% pole_specids)]

#k means on total mutation counts
kMeansMutsNoPole <- kmeans(totalTriCountsNoPole,2)
#cluster label for hypermutant (cluster with larger muts) is either 1 or 2,
#but changes per evaluation
kMeansNonHypermutLabel <- which.min(kMeansMutsNoPole$centers) 
idsNoHypermut <- names(kMeansMutsNoPole$cluster[kMeansMutsNoPole$cluster==kMeansNonHypermutLabel])
maxMutCntNoHyper <- max(totalTriCountsNoPole[idsNoHypermut]) #30148
minMutCntHyper <- min(totalTriCountsNoPole[!(names(totalTriCountsNoPole) %in% idsNoHypermut)]) #81348
#50,000 serves as cut-off

#### Construct representative non-hypermutant CRC spectra #### 
#next average the tricounts for the nonHypermut as in top of pergene_driver_trimut_analysis4.R
#then save
pcawgCrcTriCountsNonHyp <- cbind(pcawgTriCounts[,c("Mutation.type","Trinucleotide")],
                                 pcawgTriCounts[,idsNoHypermut])

pcawgCrcTriProbNonHyp <- sapply(which(colnames(pcawgCrcTriCountsNonHyp) %>% grepl("Colo",.)), 
                                function(k){
                              ( pcawgCrcTriCountsNonHyp[,k]/sum(pcawgCrcTriCountsNonHyp[,k]) ) %>% return()
                              }) 

repNonHypermutSpectra <- rowMeans(pcawgCrcTriProbNonHyp) 

##need to fix names
names(repNonHypermutSpectra ) <- sapply(1:nrow(pcawgCrcTriCountsNonHyp), function(k){
  paste0(substr(pcawgCrcTriCountsNonHyp$Trinucleotide[k],1,1), "[",
         pcawgCrcTriCountsNonHyp$Mutation.type[k],"]",
         substr(pcawgCrcTriCountsNonHyp$Trinucleotide[k],3,3)) %>% 
           return()
})
repNonHypermutSpectraDf = data.frame("MutationType" = names(repNonHypermutSpectra),
                              "SpectraVals" = as.numeric(repNonHypermutSpectra) )

# reorder to patch POLE spectra
repNonHypermutSpectraDf<-repNonHypermutSpectraDf[match(repPoleSpectraDf$MutationType, 
                                            repNonHypermutSpectraDf$MutationType), ]

repNonHypermutSpectraVec = as.numeric(repNonHypermutSpectraDf$SpectraVals) %>% as.matrix(.,col = 1)
rownames(repNonHypermutSpectraVec) = repNonHypermutSpectraDf$MutationType
colnames(repNonHypermutSpectraVec) = ("Non Hypermutant CRC Spectra")

#Fig 5B upper
plotrepNonHypermutSpectra = plot_96_profile(repNonHypermutSpectraVec,ymax=.1)+
  theme(strip.text.y = element_blank(),
        text= element_text(family = "arial",
                           size= 6)
  )


if (saveplots == T){
  pltfilename = paste(plotdir,"repNonHypermutSpectra_",
                      todaysDate,".pdf",sep = "")

  ggplot2::ggsave(filename =   pltfilename , 
                  plot =plotrepNonHypermutSpectra , 
                  device = cairo_pdf, 
                  dpi = 200, 
                  width = 7,
                  height = 2, 
                  units = "in")

}

if (save_reduceddata == T){
  filenamerepNonHypermutSpectraCsv = paste(reduced_datadir,"pcawgRepNonHypermutSpectra.csv",sep="")
  write.csv(repNonHypermutSpectraDf,
            file =filenamerepNonHypermutSpectraCsv,
            row.names =  F)
}



# MSI mutation signature --------------------------------------------------
# Use MSI signature from Brunet Guasch et al https://doi.org/10.1158/0008-5472.CAN-25-0445
# see supp table S4 of Brunet Guasch

signatures_MSI=c('SBS1','SBS5','SBS15','SBS26','SBS44','SBS57')
vector_MSI_GEL=c(0.127,0.33,0.084,.073,.31,.07)
cosmicDataDir <-   paste0(projectRoot ,
                          "data/cosmic/")
cosmicSbsVals = read.table(paste0(cosmicDataDir,"COSMIC_v3_SBS_GRCh37.txt"),
                           header = T)
repMSIvals <- sapply(1:96, function(k){
  sum(cosmicSbsVals[signatures_MSI][k,] * vector_MSI_GEL)
})

repMsiSigdf <- data.frame("MutationType" = cosmicSbsVals $Type,
                        "SpectraVals" = repMSIvals)
# reorder to patch POLE spectra
repNonHypermutSpectraDf<-repMsiSigdf[match(repPoleSpectraDf$MutationType, 
                                           repMsiSigdf $MutationType), ]

if (save_reduceddata == T){
  filenamerepMsiSpectraCsv = paste0(reduced_datadir,"repMsiSpectra.csv",sep="")
  write.csv(repMsiSigdf,
            file =filenamerepMsiSpectraCsv,
            row.names =  F)
}
