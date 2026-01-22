# Fig 6 new panels, artboard


library(magrittr)



codeDir <- paste0(projectRoot,
                  "rebuttalCodeUnclean/fig6newPanelsPlotting/")

plotdir = paste0(projectRoot,"images/krasSimPlots/")
fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveplots = F

mmToInches <- function(mm){
  return(mm/25.4)
}


# relH <- mmToInches(53.7*2)
# relW <- mmToInches(180)
# 
# quartz(width = relW,height = relH)

# mat <- matrix(
#   c(1, 2, 3,   # row 1
#     4, 5, 6),  # row 2
#   nrow = 2,
#   byrow = TRUE
# )
# 
# widths <- c(1, 12/7, 2/7)   # row 1 columns
# heights <- c(1, 1)          # row heights equal
# layout(mat, widths = widths, heights = heights)

if (savePlots == T){
  relH <- mmToInches(50) #one unit 
  relW <- mmToInches(55)
  
  filenamePltPropSel<- paste0(plotdir,
                                   fileOutCommonString,
                                   "propVsSelPlot.pdf")
  cairo_pdf(file=filenamePltPropSel ,
            width=relW,
            height=relH,
            family="Arial",
            fallback_resolution=600)
  source(paste0(codeDir,"propVsSelPlot.R"))
  dev.off()
  
  filename_abcSelScatterPlot<- paste0(plotdir,
                              fileOutCommonString,
                              "abcSelScatterPlot.pdf")
  cairo_pdf(file=filename_abcSelScatterPlot ,
            width=2*relW,
            height=relH,
            family="Arial",
            fallback_resolution=600)
  source(paste0(codeDir,"abcSelScatterPlot.R"))
  dev.off()
  
  
  filename_posteriorDeltaPlot<- paste0(plotdir,
                                      fileOutCommonString,
                                      "posteriorDeltaPlot.pdf")
  cairo_pdf(file=filename_posteriorDeltaPlot ,
            width=relW,
            height=relH,
            family="Arial",
            fallback_resolution=600)
  source(paste0(codeDir,"posteriorDeltaPlot.R"))
  dev.off()
  
  
  filename_propVsBiasPlot<- paste0(plotdir,
                                      fileOutCommonString,
                                      "propVsBiasPlot.pdf")
  cairo_pdf(file=filename_propVsBiasPlot ,
            width=relW,
            height=relH,
            family="Arial",
            fallback_resolution=600)
  source(paste0(codeDir,"propVsBiasPlot.R"))
  dev.off()
  
  
  filename_postPredMsiDens<- paste0(plotdir,
                                      fileOutCommonString,
                                      "postPredMsiDens.pdf")
  cairo_pdf(file=filename_postPredMsiDens ,
            width=relW,
            height=relH,
            family="Arial",
            fallback_resolution=600)
  source(paste0(codeDir,"postPredMsiDens.R"))
  dev.off()
  
  
  filename_postDensMuWeak<- paste0(plotdir,
                                    fileOutCommonString,
                                    "postDensMuWeak.pdf")
  cairo_pdf(file=filename_postDensMuWeak ,
            width=2*relW,
            height=relH,
            family="Arial",
            fallback_resolution=600)
  source(paste0(codeDir,"postDensMuWeak.R"))
  dev.off()
  
}



source(paste0(codeDir,"propVsSelPlot.R"))
source(paste0(codeDir,"abcSelScatterPlot.R"))
source(paste0(codeDir,"posteriorDeltaPlot.R"))
source(paste0(codeDir,"propVsBiasPlot.R"))
source(paste0(codeDir,"postPredMsiDens.R"))
# 

