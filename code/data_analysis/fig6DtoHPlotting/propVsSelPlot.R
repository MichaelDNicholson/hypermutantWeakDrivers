

library(magrittr)


mmToInches <- function(mm){
  return(mm/25.4)
}
# 
# relH <- mmToInches(53.7)
# relW <- mmToInches(180/3)
# quartz(width = relW,height = relH)

par(mar = c(4, 4, 1, 1),family = "Arial",cex = .6,las = 1)  



reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")

reduced_simdatadir <-paste(projectRoot,"reduced_data/simData/",sep = "")

fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOut = F

dfProb3rdPgridSel <- readRDS(paste0(reduced_simdatadir,"20260111.simOutWeakPropVsWeakSel.rds"))



# Relevant data parameters
biasweakNonhyp <- 8.8
biasweakPole <- 35
mutMultPole <- 100
mu3ANonhyp <- 10^(-5)


dataPweakPole <- .74 #proportion non codon 12 13 in Pole tumours
dataPweakNonHyp <- .15 #proportion non codon 12 13 muts in nonhypermutant tumours




# Plot proportions vs selection -------------------------------------------

# What is the below plot? 
# fix a muwd in nonHyp, then all rest is fixed, then fix strong driver select coef
# and find where non-hyp and 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7","purple")
strong_value1 <- 4
strong_value2<- 15
plot(subset(dfProb3rdPgridSel,ssd==4)$swd,
     subset(dfProb3rdPgridSel,ssd==4)$probsNonHyp,ylim = c(0,1.05),
     type = 'l',col =cbbPalette[4],lty = 'dashed',lwd = 2,
     ylab = "Proportion with weak driver",
     xlab ="")
lines(subset(dfProb3rdPgridSel,ssd==4)$swd,subset(dfProb3rdPgridSel,ssd==4)$probsPole,
      col = cbbPalette[9],lty = 'dashed',lwd = 2)
abline(h=dataPweakPole,col = cbbPalette[9])
abline(h= dataPweakNonHyp,col = cbbPalette[4])
lines(subset(dfProb3rdPgridSel,ssd==15)$swd,subset(dfProb3rdPgridSel,ssd==15)$probsNonHyp,ylim = c(0,1),
      type = 'l',col = cbbPalette[4],lty = 'dotted',lwd = 2)
lines(subset(dfProb3rdPgridSel,ssd==15)$swd,subset(dfProb3rdPgridSel,ssd==15)$probsPole,
      col = cbbPalette[9],lty = 'dotted',lwd = 2)
mtext( expression(paste("Weak KRAS selection coefficient (",s["weak"],')')),
      side = 1,
      line = 3,
      adj = 1,
      cex = .6) 
# Add legend showing only s[strong] = value
legend( x = -.1, y = 1.13, 
       legend = c(
         bquote(s[strong] == .(strong_value1)),
         bquote(s[strong] == .(strong_value2))
       ),
       lty = c(2, 3),  # 3 = dotted, 2 = dashed
       lwd = 2,
       bty = "n")
text(x = par("usr")[2]-.1,            # slightly inside the right edge
     y = dataPweakPole - 0.05,            # slightly below the line
     labels = expression(italic(POLE)*"-mutant"),
     col =cbbPalette[9],
     adj = c(1, 0.5))  # align text to the right of x-coordinate

text(x = par("usr")[2]-.1 , 
     y = dataPweakNonHyp - 0.05,         # slightly below the line
     labels = "Non-hypermutant",
     col =cbbPalette[4],
     adj = c(1, 0.5))
