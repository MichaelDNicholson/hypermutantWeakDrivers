# Plot characteristic of POLE mut processes
# Mut burden and signautre exposures
# Use pcawg data.

# Fig 1b,c 




# Preamble ----------------------------------------------------------------


library(magrittr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(scales)



pcawgDataDir <-  paste0(projectRoot ,
                        "data/pcawg_data/")
reduced_datadir <- paste0(projectRoot,"reduced_data/")
plotdir =  paste0(projectRoot,"images/pcawg_trinuc_plots/")
todaysDate <- gsub("-","",Sys.Date()) #for filenames



# Import/annotate data ----------------------------------------------------


#### PCAWG samples with pathogenic POLE ####
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


#### PCAWG trinuclueotide snv counts ####

#available here: https://www.synapse.org/#!Synapse:syn11726620
pcawgTriCounts =read.csv(unz(paste0(pcawgDataDir,
                                    "WGS_PCAWG_2018_02_09.zip"),
                             "WGS_PCAWG.96.csv"),header = T)
pcawgCrcTriCounts <- cbind(pcawgTriCounts[,c(1,2)],
                           pcawgTriCounts[,grepl("ColoRect",colnames(pcawgTriCounts))])

#downloaded from https://www.synapse.org/#!Synapse:syn11804040
sigExposures = read.csv(paste0(pcawgDataDir,
                               "PCAWG_sigProfiler_SBS_signatures_in_samples.csv"),
                        header = T)


##### Exclude hypermutant tumours (but non pathogenic Pole) ####
#define hypermutant tumours based on k-means clustering of total count

totalTriCounts <-  colSums(pcawgCrcTriCounts[,3:ncol(pcawgCrcTriCounts)])
TriCountsSpec <- names(totalTriCounts) %>% gsub("ColoRect.AdenoCA..","",.)

#exclude already identified pole edm containing crcs
totalTriCountsNoPole <- totalTriCounts[!(TriCountsSpec %in% pole_specids)]

# Classify as hypermutant or not
#k means on total mutation counts
kMeansMutsNoPole <- kmeans(totalTriCountsNoPole,2)
kMeansNonHypermutLabel <- which.min(kMeansMutsNoPole$centers) 
idsNoHypermut <- names(kMeansMutsNoPole$cluster[kMeansMutsNoPole$cluster==kMeansNonHypermutLabel])




# Mutation burden plots ------------------------------------------------
totalTriCountsWPole <- totalTriCounts[(TriCountsSpec %in% pole_specids)]
totalTriCountsNoHyp <- totalTriCounts[ idsNoHypermut]

#10^6/(3*10^9) is approximately per sbs per megabase
MutsPerMegaNoHyp <- (totalTriCountsNoHyp/(3*10^3)) 
MutsPerMegaWPole <- (totalTriCountsWPole/(3*10^3)) 

MutsPerMegaNoHypDf <- data.frame("Mut" =sort(MutsPerMegaNoHyp),
                                 "index"= 1:length(MutsPerMegaNoHyp))
MutsPerMegaWPoleDf <- data.frame("Mut" =sort(MutsPerMegaWPole),
                                 "index"= 1:length(MutsPerMegaWPole))


#### make plot just to get legend ####
MutsPerMegaBoth <- data.frame("Mut" = sort(c(MutsPerMegaNoHyp,
                                             MutsPerMegaWPole)),
                              "index" = 1:length(c(MutsPerMegaNoHyp,
                                                   MutsPerMegaWPole)),
                              "type" = c(rep("Non-hypermutant",length(MutsPerMegaNoHyp)),
                                         rep("With POLE\nmutation",length( MutsPerMegaWPole)))) 

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

pairedCols <- brewer.pal(3, "Paired")

pltForLegend <- ggplot(MutsPerMegaBoth,aes(x=index,y=Mut,color = type))+
  geom_point()+
  scale_color_manual(values =pairedCols[c(2,1)],name = "",
                     labels= c("Non-hypermutant",
                               expression(italic("POLE")~"mutated")))
goodLegend <- get_legend(pltForLegend)

#### burden by non hypermutant and those with pole ####
pltNonHyp <- ggplot(MutsPerMegaNoHypDf, aes(x= index,y=Mut))+
  geom_point(size = .25,color = pairedCols[2])+
  geom_hline(yintercept = mean(MutsPerMegaNoHyp),
             linetype = "dashed")+
  scale_y_log10(breaks = c(1,10,100,1000), 
                limits = c(1,1000))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(family = "arial",
                            size= 6))+
  ylab("Single base substitutions\nper megabase")

pltWPole <- ggplot(MutsPerMegaWPoleDf, aes(x= index,y=Mut))+
  geom_point(size=.25,color = pairedCols[1])+
  geom_hline(yintercept = mean(MutsPerMegaWPole),
             linetype = "dashed")+
  scale_y_log10(breaks = c(1,10,100,1000), 
                limits = c(1,1000))+
  theme_classic()+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line.y = element_blank(),
        text = element_text(family = "arial",
                            size= 6))
  

plotBoth <- plot_grid(pltNonHyp,pltWPole ,
                      rel_widths = c(nrow(MutsPerMegaNoHypDf),
                                     2*nrow(MutsPerMegaWPoleDf)) )

plotBothLegend <- plot_grid(plotBoth,goodLegend,rel_widths = c(1,1))



if (saveplots == T){
  pltfilename = paste(plotdir,"mutBurdenNonHypPoleNoLeg_",
                      todaysDate,".pdf",sep = "")

  ggplot2::ggsave(filename =   pltfilename , 
                  plot =  plotBoth , 
                  device = cairo_pdf, 
                  dpi = 300, 
                  width = 50,
                  height = 45, 
                  units = "mm")
  

}


# Exposure vectors --------------------------------------------------------

# Average sig values for POLE (av exposure >0.01)
sigExposuresPoleSamples = subset(sigExposures  ,Sample.Names %in% 
                                    pole_specids)
indexSbsCols <- colnames(sigExposuresPoleSamples ) %>% grepl("SBS",.) %>% which()
sigExposureProbsPole = lapply(1:nrow(sigExposuresPoleSamples), function(k){
  return(sigExposuresPoleSamples[k,indexSbsCols]/
           sum(sigExposuresPoleSamples[k,indexSbsCols]))
}) %>% do.call(rbind,.)
sigExposureProbsPoleMean <- colMeans(sigExposureProbsPole)
sigExposureProbsPoleMean01 <- sigExposureProbsPoleMean[sigExposureProbsPoleMean>.01]

# Average sig values for non-hypermutant (av exposure  > 0.01)

sampIdsNoHypermut <- sapply(idsNoHypermut, function(x){
  gsub(pattern = "ColoRect.AdenoCA..","",x) %>% return()  
})
sigExposuresNoPoleSamples <- subset(sigExposures  ,Sample.Names %in% 
                                      sampIdsNoHypermut)
sigExposureProbsNoPole = lapply(1:nrow( sigExposuresNoPoleSamples), function(k){
  return(sigExposuresNoPoleSamples[k,indexSbsCols]/
           sum(sigExposuresNoPoleSamples[k,indexSbsCols]))
}) %>% do.call(rbind,.)
sigExposureProbsNoPoleMean <- colMeans(sigExposureProbsNoPole)
sigExposureProbsNoPoleMean01 <-sigExposureProbsNoPoleMean[sigExposureProbsNoPoleMean>.01]


# Plot av sig exposure for Pole/non-pole samples
commonSigsSeen01 <- union(names(sigExposureProbsNoPoleMean[sigExposureProbsNoPoleMean>0.01]),
                          names(sigExposureProbsPoleMean[sigExposureProbsPoleMean>0.01])  )
commonSigsSeen01Ordered <- c("SBS1",
                             "SBS5", 
                             "SBS10a",
                             "SBS10b",
                             "SBS17b",
                             "SBS18",
                             "SBS28",
                             "SBS37",
                             "SBS40",
                             "SBS44")

statusFactor <- factor(c(rep("POLE-mutant",length(commonSigsSeen01Ordered)),
                         rep("Non-hypermutant",length(commonSigsSeen01Ordered))),
                       levels = c("POLE-mutant","Non-hypermutant"))
exposureNormDf <- data.frame("sig" = c(commonSigsSeen01Ordered,commonSigsSeen01Ordered),
                             "normExpos" = c(sigExposureProbsPoleMean[commonSigsSeen01Ordered],
                                             sigExposureProbsNoPoleMean[commonSigsSeen01Ordered]),
                             "status" = statusFactor )
exposureNormDf$sigFact <- factor(exposureNormDf$sig,
                                 level= commonSigsSeen01Ordered)
plotExposBars <- ggplot(data=exposureNormDf, aes(x=sigFact, y=normExpos, fill=status)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title=element_blank(),
        legend.position = "top")+
  ylab("Average signature proportion\n(PCAWG CRCs)")+
  xlab("Mutational signature")





if (saveplots == T){
  
  mmToInches = function(mm){return(mm*0.0393701)}
  hwratio = 2.5/4
  
  pltfilename = paste(plotdir,"mutExposuresNonHypPole_",
                      todaysDate,".pdf",sep = "")
  save_plot(filename = pltfilename,
            plotExposBars ,
            base_height = 2.8,
            base_width = 3.5)
}

