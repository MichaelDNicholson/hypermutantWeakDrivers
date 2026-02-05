
library(magrittr)


mmToInches <- function(mm){
  return(mm/25.4)
}

# relH <- mmToInches(53.7)
# relW <- mmToInches(180/3*2)
# quartz(width = relW,height = relH)




reduced_datadir <-paste(projectRoot,"reduced_data/",sep = "")

reduced_simdatadir <-paste(projectRoot,"reduced_data/simData/",sep = "")

fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")

saveOut = F


#simulation functions
source(paste0(projectRoot,"code/modelling/simulationFunctionsClean.R"))


# Relevant data parameters
biasweakNonhyp <- 8.8
biasweakPole <- 35
mutMultPole <- 100
mu3ANonhyp <- 10^(-5)


dataPweakPole <- .74 #proportion non codon 12 13 in Pole tumours
dataPweakNonHyp <- .15 #proportion non codon 12 13 muts in nonhypermutant tumours



# Thresholds from eq !
deltaThreshPole <- -1*log10(biasweakPole)/log10(mu3ANonhyp*mutMultPole)
deltaThreshNonHyp <- -1*log10(biasweakNonhyp)/log10(mu3ANonhyp)






# ABC accepted selection params  -------------------------------------------

dfParamsProps <- readRDS(paste0(reduced_simdatadir,"20260108.simKrasInf50Ksims.rds"))


abcEps <- .05
dfParamsPropsAbcPass <- subset(dfParamsProps,distTotal<abcEps)

print(paste0("mean fold change increase growth rate ",
       as.character(mean( (1+dfParamsPropsAbcPass $ssd)/(1+dfParamsPropsAbcPass$swd)))))
     
print(paste0("CI95 ",
             as.character(quantile( (1+dfParamsPropsAbcPass $ssd)/(1+dfParamsPropsAbcPass$swd),
                                    c(0.025,.975)))))



# scatter plot of accepted selection params -------------------------------

x <- dfParamsProps$swd
y <- dfParamsProps$ssd
z <- dfParamsProps$delta
d <- dfParamsProps$distTotal



# Keep only complete cases
ok <- complete.cases(x, y, z, d)
x <- x[ok]; y <- y[ok]; z <- z[ok]; d <- d[ok]

# -----------------------------
# Colour palette
# -----------------------------
ncol <- 100
cols <- colorRampPalette(c("purple4", "orchid", "gold", "darkorange"))(100)

# Map delta to colours
z_cut <- cut(z, breaks = ncol, labels = FALSE)
col_vec <- cols[z_cut]

# Shape based on distTotal
pch_vec <- ifelse(d >= abcEps, 1, 16)

# Reduce opacity for open circles
alpha <- 0.05
col_vec_alpha <- col_vec
col_vec_alpha[pch_vec == 1] <- rgb(t(col2rgb(col_vec[pch_vec == 1])/255), alpha = alpha)

# -----------------------------
# Layout: scatter + color bar
# -----------------------------
layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
# par(fig = c(0, 6/7, 0, 1), new = TRUE)
# -----------------------------
# Scatter plot
# -----------------------------
# oldpar <- par(no.readonly = TRUE)
# par(fig = c(0, 6/7, 0, 1))

par(mar = c(4, 4.4, 2, 1),family = "Arial", cex = 0.6,
    las = 1)


plot(
  x, y,
  xlim = c(0, 3),
  ylim = c(0, 20),
  col = col_vec_alpha,
  pch = pch_vec,
  xlab = expression(paste("Weak selection coefficient (",s["weak"],')')),
  ylab = expression(paste("Strong selection coefficient (",s["strong"],')'))
)

# Shape legend
legend(
  "bottomright",
  legend = c("ABC accepted samples", "ABC rejected samples"),
  pch = c(16, 1),
  col = c("black", "black"),
  bty = "n"
)

# -----------------------------
# Threshold curves
# -----------------------------
swd_seq <- seq(min(x), max(x), length.out = 500)

# Function to compute ssd from swd and delta
ssd_from_delta <- function(swd, delta) {
  ssd <- (swd + delta*(1 + swd)) / (1 - delta*(1 + swd))
  # Only keep ssd > swd
  ssd[ssd <= swd] <- NA
  return(ssd)
}

ssd_threshNonHyp <- ssd_from_delta(swd_seq, deltaThreshNonHyp)
ssd_threshPole   <- ssd_from_delta(swd_seq, deltaThreshPole)
ssd_threshMeanDelta   <- ssd_from_delta(swd_seq, deltaMean)


# Overlay curves (omit NA automatically)
lines(swd_seq, ssd_threshNonHyp, col = "black", lwd = 2, lty = 2)   # dotted
lines(swd_seq, ssd_threshPole, col = "black", lwd = 2, lty = 2)       # dashed
lines(swd_seq,ssd_threshMeanDelta, col = "black", lwd = 2, lty = 1)       # dashed


# -----------------------------
# Color bar legend
# -----------------------------
# -----------------------------
# Color bar legend
# -----------------------------
# par(fig = c(6/7, 1, 0, 1), new = TRUE)

par(mar = c(2, 1, 2, 4),family = "Arial", cex = 0.6,
    las = 0)

z_seq <- seq(min(z), max(z), length.out = ncol)

image(
  x = c(0, 1),
  y = z_seq,
  z = matrix(z_seq[-length(z_seq)], nrow = 1),
  col = cols,
  axes = FALSE,
  xlab = "",
  ylab = ""
)

axis(4)
box()

# Colour bar label with symbol
mtext(
  expression("Scaled selection difference (" * Delta * ")"),
  side = 4,
  line = 2.5,
  cex = 0.6,
)
# par(oldpar)


# dev.off()
# graphics.off()
