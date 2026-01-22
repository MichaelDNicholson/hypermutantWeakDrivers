# Script to determine mutational bias to individual sites
# and groups of drivers under POLE/non-hypermutant CRC
# mutational processes

# Preamble ----------------------------------------------------------------
library(magrittr)
library(ggplot2)
library(cowplot)


reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")
driverListDir <- paste(projectRoot,"data/downloaded_gene_data/",sep = "")
plotDir <- paste(projectRoot,"images/krasMutationalBias/",sep = "")
fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOutput <- F

##for plotting
clone_colors = c("#619CFF",
                  "darkgoldenrod1",
                  "#F8766D" )
                  


#load in dataframe of kras mutations, annotated with driver info
combinedKrasList <- read.csv(file =paste0(reduced_datadir,"combinedKrasListDriverAnnotated.csv"),
                             header = T)
combinedKrasList$mutBiasPoleNonPole <- combinedKrasList$relMutRateCrcPole/
                                      combinedKrasList$relMutRateCrc

combinedKrasList$mutAACDS <- sapply(1:nrow(combinedKrasList), function(k){
  paste(combinedKrasList$mutAA[k],
        combinedKrasList$mutCds[k],
        sep = ", ") %>% return
})




combinedKrasList$rowIndex <- 1:nrow(combinedKrasList)


# Functions ---------------------------------------------------------------



getRelProbWeakStrong <- function(df,class){
#df is combinedKrasList 
#function returns relative prob of weak vs strong driver under 
#pole and non-pole processes
  
  pweakPOLE <- sum(df$relMutRateCrcPole[df[,class]=="Weak"])
  pstrongPOLE <- sum(df$relMutRateCrcPole[df[,class]=="Strong"])
  
  pweakNonHyp <- sum(df$relMutRateCrc[df[,class]=="Weak"])
  pstrongNonHyp <- sum(df$relMutRateCrc[df[,class]=="Strong"])
  
  return(c(pweakPOLE/pstrongPOLE,pweakNonHyp/pstrongNonHyp))
}


  
  getAllPlotsGivenClass <- function(classVar){
  #for a given way of decidign strong vs
  #weak kras driver (classVar), this function returns 
  #all relevant plots
          
        plotFoldchange = 
          ggplot(data=combinedKrasList ,aes(color = !!sym(classVar))) +
          xlab("")+
          ylab("Ratio mutation proportions \n (POLE/non-hypermutant)")+
          geom_segment(aes(y=0,yend=relMutRateCrcPole/relMutRateCrc,
                           x=rowIndex, xend=rowIndex),
                       size = 3)+
          theme_classic()+
          theme( legend.position = c(.1,.9),
                 legend.title = element_blank(),
                 legend.text = element_text(size=12),
                 axis.text.x = element_text(angle = 90, 
                                            vjust = 0.5, 
                                            hjust=1)) +
          scale_color_manual(values = clone_colors[c(3,2)])+
          scale_x_continuous(labels = as.character(combinedKrasList$mutAACDS ),
                             breaks = 1:nrow(combinedKrasList ))+
          scale_y_continuous(breaks = 0:5)+
          geom_hline(yintercept = 1, linetype = "dashed",linewidth = 1)
       
        
        
        relRates<- getRelProbWeakStrong(combinedKrasList,classVar)
        probRatios= data.frame("relprob"= relRates ,
                                   "Background" = c("POLE-mutant","non-hypermutant"),
                                   "pos" = c(1,2))
        
        pltRatioRates <- ggplot(probRatios , aes(x=Background, y=relprob)) + 
          geom_bar(stat = "identity",
                   color="black", 
                   fill=rgb(0.1,0.4,0.5,0.7),alpha=1,
                   width = .5)+
          theme_classic()+
          xlab("Mutational process")+
          ylab(expression(frac("KRAS weak-driver mutation rate",
                               "KRAS strong-driver mutation rate")))+
          geom_text(aes(label=signif(relprob,2)),
                    vjust=1.6, color="white", size=3.5)
        
        
        

        
        
      
        
        return(list(plotFoldchange ,
                    pltRatioRates))
  }
  
  
  # Make plots and get mutation bias ratios --------------------------------------------------------------
  
  
  
  possibleClassVars <- colnames(combinedKrasList)[colnames(combinedKrasList) %>%
                                                    grep("class",.)]
  mutBiasPoleNonPoleByClass <- sapply(possibleClassVars, function(classVar){
    relRates<- getRelProbWeakStrong(combinedKrasList,classVar)
    return(c(relRates[1]/relRates[2]))
  })
  #Ratio of bias weak over strong for Pole/non-Pole
  mutBiasPoleNonPoleByClassDf <- data.frame("RatioBiasPoleNonPole" =mutBiasPoleNonPoleByClass,
                                            "ClassCriteria" = names(mutBiasPoleNonPoleByClass))
  
  
  possibleClassVarsExplanation <- c(
    "Kmeans on bioactivation score from Haigis 17",
    "Mutations in codons 12, 13 are strong drivers",
    "Kmeans on observed:expected variant counts in cBioPortal",
    "Kmeans on observed:expected variant counts in cBioPortal annotated as MSS"
    
  )
  names(possibleClassVarsExplanation) <- possibleClassVars
  
  listRelPlotsByClassVar <- lapply(possibleClassVars,getAllPlotsGivenClass )
  
  listRelPlotsByClassVarGrid <- lapply(1:length(possibleClassVars), function(k){
    classTitle <- ggdraw() + 
      draw_label(paste0("Assigning driver class by: ",
                  possibleClassVarsExplanation[k]),
                 fontface='bold')
    gridRepPlots <- plot_grid(plotlist = listRelPlotsByClassVar[[k]],
                              nrow = 1)
    plot_grid(classTitle,gridRepPlots,
              rel_heights =c(.2,1),ncol = 1) %>% return
    
  })
  
  allRelPlotsByClass <- plot_grid(plotlist  = listRelPlotsByClassVarGrid ,
                                  nrow = length(listRelPlotsByClassVar))
  
  #Fig 5C,D
  exon2VariantFoldChange <- listRelPlotsByClassVar[[2]][[1]]+ggtitle("")
  exon2MutBiasFoldChange <- listRelPlotsByClassVar[[2]][[2]]+ggtitle("")
  

# Save plots and mutational biases --------------------------------------------------------------

if (saveOutput == T){
  pltfilename <- paste0(plotDir,"allMutBiasPlotsByClass.pdf")
  
  save_plot(pltfilename,allRelPlotsByClass,
            base_height = 20, base_asp = 1.3)
  
  pltfilename <- paste0(plotDir,"cBioCntsVariantFoldChange.pdf")
  save_plot(pltfilename,cBioCntsVariantFoldChange ,
            base_height = 3.5, base_asp = 2)
 

  
  pltfilename <- paste0(plotDir,fileOutCommonString,"exon2VariantFoldChange.pdf")
  save_plot(pltfilename,exon2VariantFoldChange ,
            base_height = 3.5, base_asp = 2.4)
  

  
  filenameMutBias <- paste0(reduced_datadir,"ratioMutBiasPoleNonHyp.csv")
  
  write.csv(mutBiasPoleNonPoleByClassDf,
            file = filenameMutBias,
            row.names =  F)
  
  
  
}

  
  

# MSI section ---------------------------------------------------------
  # combinedKrasList <- read.csv(file =paste0(reduced_datadir,"combinedKrasListDriverAnnotated.csv"),
  #                              header = T)
  combinedKrasList$mutBiasMSINonPole <- combinedKrasList$relMutRateMSI/combinedKrasList$relMutRateCrc
    
  
  
  getRelProbWeakStrongMSI <- function(df,class){
    #df is combinedKrasList 
    #function returns relative prob of weak vs strong driver under 
    #pole and non-pole processes
    
    pweakMSI <- sum(df$relMutRateMSI[df[,class]=="Weak"])
    pstrongMSI <- sum(df$relMutRateMSI[df[,class]=="Strong"])
    
    pweakNonHyp <- sum(df$relMutRateCrc[df[,class]=="Weak"])
    pstrongNonHyp <- sum(df$relMutRateCrc[df[,class]=="Strong"])
    
    return(c(   pweakMSI/ pstrongMSI,pweakNonHyp/pstrongNonHyp))
  }
  
  mutBiasMsiNoHypByClass <- sapply(possibleClassVars, function(classVar){
    relRates<-   getRelProbWeakStrongMSI(combinedKrasList,classVar)
    return(c(relRates[1]/relRates[2]))
  })
  
  mutBiasMSINonPoleByClassDf <- data.frame("RatioBiasMSINonPole" =  mutBiasMsiNoHypByClass,
                                            "ClassCriteria" = names(  mutBiasMsiNoHypByClass))
  
  
  
  
  classVar <-  "class1213"  
  ggplot(data=combinedKrasList ,aes(color = !!sym(classVar))) +
    xlab("")+
    ylab("Ratio mutation proportions \n (POLE/non-hypermutant)")+
    geom_segment(aes(y=0,yend=mutBiasMSINonPole,
                     x=rowIndex, xend=rowIndex),
                 size = 3)+
    theme_classic()+
    theme( legend.position = c(.1,.9),
           legend.title = element_blank(),
           legend.text = element_text(size=12),
           axis.text.x = element_text(angle = 90, 
                                      vjust = 0.5, 
                                      hjust=1)) +
    scale_color_manual(values = clone_colors[c(3,2)])+
    scale_x_continuous(labels = as.character(combinedKrasList$mutAACDS ),
                       breaks = 1:nrow(combinedKrasList ))+
    scale_y_continuous(breaks = 0:5)+
    geom_hline(yintercept = 1, linetype = "dashed",linewidth = 1)
  
  
