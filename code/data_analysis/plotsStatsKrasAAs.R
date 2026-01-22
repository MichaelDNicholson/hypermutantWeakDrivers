# Take in KRAS mutations annotated by POLE status and 
# extract statistics/make plots




library(magrittr)
library(data.table)
library(plyr)
library(ggplot2)


# Preamble ----------------------------------------------------------------


reduced_datadir <- paste0(projectRoot,"reduced_data/")
plotdir = paste0(projectRoot,"images/krasPublicData/")

fileOutCommonString <- paste0(format(Sys.time(),"%Y%m%d"),".")
saveplots = F

# load in kras mutatation annotated.
krasDataMisWPoleEdm <- read.csv( file = paste0(reduced_datadir,"cbioKrasDataWPoleEdm.csv"))
krasDataMisNoPoleEdm <- read.csv( file = paste0(reduced_datadir,"cbioKrasDataNoPoleEdm.csv"))


#check if samples are duplicates
#only keep patient  mutation combos (some patients have multiple samples)
noDupsKrasMisWPoleEdm <- ddply(krasDataMisWPoleEdm[,c("PatientID","Protein.Change","aaNumber")],
                               .(PatientID,`Protein.Change`,aaNumber),nrow)
noDupsKrasMisNoPoleEdm <- ddply(krasDataMisNoPoleEdm[,c("PatientID","Protein.Change","aaNumber")],
                                .(PatientID,`Protein.Change`,aaNumber),nrow)


# Proportion of kras codons mutated overall -------------------------------
#for statements like codon 12 mutated in x% of crcs.
# Figure 1a
nonpoleedm_kras_aa = noDupsKrasMisNoPoleEdm$aaNumber %>% as.numeric
poleedm_kras_aa =noDupsKrasMisWPoleEdm$aaNumber %>% as.numeric
both_kras_aa = c(nonpoleedm_kras_aa,poleedm_kras_aa)
nmuts_pole_kras = length(poleedm_kras_aa)


codonMutProbNonPole <- table(nonpoleedm_kras_aa)/length(nonpoleedm_kras_aa)
codonMutProbPole <- table(poleedm_kras_aa)/length(poleedm_kras_aa)


codonMutPropAll <- (both_kras_aa %>% table)/length(both_kras_aa )

prop1213NoPole <- sum(nonpoleedm_kras_aa %in% c(12,13))/length(nonpoleedm_kras_aa)
prop1213Pole <- sum(poleedm_kras_aa %in% c(12,13))/length(poleedm_kras_aa)

dfcodon1213 <- data.frame(prop = c(prop1213NoPole,1-prop1213NoPole ,
                                   prop1213Pole,1-prop1213Pole),
                          type = c("No POLE-edm","No POLE-edm",
                                   "With POLE-edm","With POLE-edm" ),
                          codons = c("Codon 12/13", "Outside codon 12/13",
                                     "Codon 12/13", "Outside codon 12/13"))

dfcodon1213$typeFactor <- factor(dfcodon1213$type,levels= c("With POLE-edm",
                                                            "No POLE-edm"))

plotCodon1213Cbio <- ggplot(data=dfcodon1213 , aes(x=codons, y=prop, fill=typeFactor)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position = 'top',
        legend.justification='left',
        text = element_text(family = "arial",
                            size = 6))+
  ylab(expression(atop(paste("Proportion of ",italic("KRAS"))," variants in CRC")))+
  xlab("")


#test: association of POLE tumours and non-codon 12/13 KRAs
fisher.test(rbind(c(sum(nonpoleedm_kras_aa %in% c(12,13)),
                    sum(!(nonpoleedm_kras_aa %in% c(12,13)))),
                  c(sum(poleedm_kras_aa %in% c(12,13)),
                    sum(!(poleedm_kras_aa %in% c(12,13))))))


if (saveplots == T){
  pltfilename = paste(plotdir,fileOutCommonString,
                      "plotCodon1213Cbio.pdf",sep = "")
  whratio = .9
  w=52
  ggplot2::ggsave(filename =   pltfilename , 
                  plot =  plotCodon1213Cbio , 
                  device = cairo_pdf, 
                  dpi = 200, 
                  width = w,
                  height = w/whratio, 
                  units = "mm")
  system(paste0('open "', pltfilename, '"'))
  # save_plot(filename= pltfilename,plotCodon1213Cbio,
  #           base_height = 3,base_asp = 1.2)
  
}

# Test if mutated residue in kras is sig. different in pole/non-pole samples --------------------------
#permutation test on mean residue
# Figure 5A



nruns = 10^5
empirical_dist_mean_aaval  = sapply(1:nruns, function(k){
  return(mean(sample(both_kras_aa ,nmuts_pole_kras,replace = F )))
  
}) %>% as.data.frame()

pval = max(1/nruns,sum(empirical_dist_mean_aaval>= mean(poleedm_kras_aa))/nruns)


ggplot(empirical_dist_mean_aaval , aes(x=.))+
  geom_histogram(color="darkblue", fill="lightblue",binwidth = 2)+
  theme_bw()+
  geom_vline(xintercept = mean(poleedm_kras_aa),color = "red")+
  xlab("Mean residue number per subsample")+
  ylab("Count")


aa_val_histcnt_nonpole = sapply(1:188,function(k){
  sum(nonpoleedm_kras_aa == k) %>% return()
})
aa_val_histcnt_pole = sapply(1:188,function(k){
  sum(poleedm_kras_aa  == k) %>% return()
})


df_aavalhist = data.frame("aa_val" = 1:188,
                          "count_pole" =aa_val_histcnt_pole ,
                          "count_nonpole" = aa_val_histcnt_nonpole,
                          "prop_pole" = aa_val_histcnt_pole/length(poleedm_kras_aa ),
                          "prop_nonpole" = aa_val_histcnt_nonpole/length(nonpoleedm_kras_aa ))


codons_seen_nonpole = nonpoleedm_kras_aa %>% unique
codons_seen_pole = poleedm_kras_aa %>% unique
codons_seeneither = union(codons_seen_nonpole,codons_seen_pole) %>% sort
prop_aaval_atcodonseen_pole = sapply(codons_seeneither, function(x){
  subset(df_aavalhist,aa_val == x)$prop_pole %>% return
})
prop_aaval_atcodonseen_nonpole = sapply(codons_seeneither, function(x){
  subset(df_aavalhist,aa_val == x)$prop_nonpole %>% return
})


df_aaval_propboth = data.frame("aa_val" = c(codons_seeneither,codons_seeneither) ,
                               "prop" = c(prop_aaval_atcodonseen_nonpole,
                                          prop_aaval_atcodonseen_pole ),
                               "status" = c(rep("POLE wild type",length(codons_seeneither)),
                                            rep("POLE-mutant",length(codons_seeneither))))

df_aaval_propboth$aa_val_factor = factor(df_aaval_propboth$aa_val %>% as.character,
                                         levels=unique(codons_seeneither %>% as.character()))

df_aaval_propboth$g01inEither <- sapply(df_aaval_propboth$aa_val, function(x){
  relsub <- subset(df_aaval_propboth, aa_val == x)
  aaG01inOnetype <- sum(relsub$prop>.01)
  return(!(aaG01inOnetype==0))
})

annotate_str = paste("POLE wild type n = ",length(nonpoleedm_kras_aa) %>% as.character(),
                     "\n",
                     "POLE-mutant n = ",length(poleedm_kras_aa) %>% as.character(),
                     "\n",
                     "p<", as.character(pval ),", permutation test",sep = "")


plot_aaval_propboth_g01 = ggplot(subset(df_aaval_propboth,
                                        g01inEither == T),
                                 aes(x=aa_val_factor, y=prop, fill=status)) +
  geom_bar(stat="identity", position=position_dodge(),
           color = "black")+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  xlab(expression(paste(italic("KRAS")," mutated residue")))+
  ylab("Proportion of mutations")+
  annotate(geom="text", x=length(codons_seeneither)-21,
           y=.5, label=annotate_str,
           hjust = 0 ,
           size = 6/.pt)+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 0),
        text = element_text(family = "arial",size = 6))

if (saveplots == T){
  pltfilename = paste(plotdir,
                      fileOutCommonString,
                      "cbioportal_aaval_propboth_g01.pdf",sep = "")
  
  ggplot2::ggsave(filename =   pltfilename , 
                  plot =plot_aaval_propboth_g01  , 
                  device = cairo_pdf, 
                  dpi = 200, 
                  width = 65,
                  height = 60, 
                  units = "mm")
  # ggplot2::ggsave(filename =   pltfilename , 
  #                 plot =plot_aaval_propboth_g01  , 
  #                 device = cairo_pdf, 
  #                 dpi = 200, 
  #                 width = 170/2,
  #                 height = 54, 
  #                 units = "mm")
  
  
}



