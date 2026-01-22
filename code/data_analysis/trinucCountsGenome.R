library(Biostrings)
library("BSgenome.Hsapiens.UCSC.hg38")


reduced_datadir = paste(projectRoot,"reduced_data/",sep = "")

saveoutput = F

#get tri counts over primary assembly (https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19)
autosome_chr = sapply(seq(1,22),function(k)paste("chr",as.character(k),sep=""))
sex_chr = c("chrX","chrY")
mit_chr = "chrM"
rands = names(Hsapiens) [names(Hsapiens) %>% grepl("random",.)]
unpl = names(Hsapiens) [names(Hsapiens) %>% grepl("chrUn",.)]
relhsap_names = c(autosome_chr,sex_chr,mit_chr,rands,unpl)
tricnts_byname = lapply(relhsap_names,function(n){
  print(n)
  return(trinucleotideFrequency(Hsapiens[[n]]))
})

tricnts_byname_mat = do.call(rbind,tricnts_byname)
tricnts_primaryassemb = colSums(tricnts_byname_mat) 
tricnts_primaryassemb_df <- data.frame("Context" = names(tricnts_primaryassemb), 
                             "Count" = as.numeric(tricnts_primaryassemb))

if (saveoutput == T){
  filename = paste0(reduced_datadir,"tricntsHg38PriAssembly.csv")
  write.csv(tricnts_primaryassemb_df ,
            file =  filename,
            row.names =  F)
}
