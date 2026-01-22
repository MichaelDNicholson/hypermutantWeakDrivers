# Repository for Hyper-mutational Processes Manuscript

This repository contains all code and processed data required to reproduce the figures and tables for the manuscript:

**Hyper-mutational processes provide a head-start for weak cancer drivers: explaining atypical KRAS variants**

To run:
- Clone repository
- Run code/setupDirectory.R
  

---

## Repository Structure

- `data/` — Raw public and reference datasets
- `reduced_data/` — Processed and intermediate data files
- `code/data_analysis/` — Statistical analyses and figure plotting
- `code/modelling/` — Analytical and simulation-based analysis


---

## Figures and Panels

### Figure 1. Atypical KRAS variants and mutational processes in POLE-mutant colorectal cancer

- **Panel A:** KRAS codon 12/13 mutations in POLE vs non-POLE CRCs  
  **Scripts:**
  - `data_analysis/annotatePublicMutData.R` — data cleaning and annotation
  - `data_analysis/plotsStatsKrasAAs.R` — statistical analysis and plotting  
  **Input data:**
  - `data/cbioportal/cbioportal_02052023_CRC_combined_study_clinical_data.tsv`
  - `data/cbioportal/cbioportal_02052023_CRC_KRAS_mutdata.tsv`
  - `data/cbioportal/cbioportal_02052023_CRC_POLE_mutdata.tsv`
  - `reduced_data/pathogenicPoleList.csv`  
  **Intermediate data:**
  - `reduced_data/cbioKrasDataWPoleEdm.csv`
  - `reduced_data/cbioKrasDataNoPoleEdm.csv`

- **Panel B:** Mutation burden and mutational signatures in POLE CRCs  
  **Script:** `data_analysis/poleMutCharacteristics.R`


### Figure 3. Mutational biases can lead to dominance of the weak driver clone

- **Panels A–C:** Clone size trajectory plots  
  **Script:** `code/modelling/doApproxSimTrajectories.R` (Clone size trajectory plots)

### Figure 4. Parameter regimes determining weak vs strong driver dominance

- **Panel A:** Probability of weak/strong clone dominance at fixed population size  
  **Script:** `code/modelling/doApproxSimsVsAnalytic.R` (Contour plots)

- **Panel B:** Difference in dominance probability versus selection coefficient  
  **Script:** `code/modelling/paramRegimesAnalytic.R` (Delta function vs selection parameters)

- **Panel C:** Regimes for emergence of a third driver via weak vs strong driver path  
  **Script:** `code/modelling/paramRegimesAnalytic.R` (Switch driver between subtypes)

### Figure 5. Altered mutational processes and non-canonical KRAS variants in POLE-mutant CRCs

- **Panel A:** Codons mutated in KRAS (POLE vs non-POLE)  
  **Script:** `data_analysis/plotsStatsKrasAAs.R`  
  **Input data:** 
  - `reduced_data/cbioKrasDataWPoleEdm.csv`
  - `reduced_data/cbioKrasDataNoPoleEdm.csv`

- **Panel B:** Mutational spectra  
  **Script:** `data_analysis/representativeSpectraMutBurden.R`  
  **Input data:** 
  - `data/pcawg_data/icgc_pcawg_crc_pole_sept21_simple_somatic_mutation.open.tsv.gz`
  - `data/pcawg_data/PCAWG_sigProfiler_SBS_signatures_in_samples.csv`
  - `data/pcawg_data/WGS_PCAWG_2018_02_09.zip`
  - `data/cosmic/COSMIC_v3_SBS_GRCh37.txt`
  - `reduced_data/pathogenicPoleList.csv`

- **Panels C–D:** Classification of KRAS variants as weak or strong drivers and mutational bias  
  **Scripts:**
  - `data_analysis/trinucCountsGenome.R`
  - `data_analysis/curateKrasDriverList.R`
  - `data_analysis/annotateKrasDriverClass.R`
  - `data_analysis/mutationalBiasAnnotation.R`  
  **Input data:**
  - `data/miscellaneous/V94_38_MUTANTCENSUS_KRAS.csv`
  - `data/miscellaneous/IntOGen-Distribution-KRAS-COREAD.tsv`
  - `data/miscellaneous/custom_KRAS_driverlist.xlsx`
  - `data/miscellaneous/Haigis18Activationscore.xlsx`
  - `reduced_data/pcawgRepPoleSpectra.csv`
  - `reduced_data/pcawgRepNonHypermutSpectra.csv`
  - `reduced_data/tricntsHg38PriAssembly.csv`
  - `reduced_data/combinedKrasList.csv`

### Figure 6. Conditions under which mutational processes alter detected KRAS drivers

- **Panels A–C:** Selection versus KRAS allele flipping between subtypes  
  **Script:** `code/modelling/paramRegimesAnalytic.R`

- **Panels D–H:** Simulation-based inference and MSI enrichment  
  **Simulation scripts:**
  - `code/modelling/simsKras.R`
  - `code/modelling/generateAbcSimData.R`  
  **Plotting scripts:**
  - `data_analysis/fig6DtoHPlotting/propVsSelPlot.R`
  - `data_analysis/fig6DtoHPlotting/abcSelScatterPlot.R`
  - `data_analysis/fig6DtoHPlotting/posteriorDeltaPlot.R`
  - `data_analysis/fig6DtoHPlotting/propVsBiasPlot.R`
  - `data_analysis/fig6DtoHPlotting/postPredMsiDens.R`
  - `data_analysis/enrichmentWeakStrongMSI.R`  
  **Input data:**
  - `reduced_data/contigencyTables.rds`

---

## Supplementary Figures and Tables

- **Supplementary Figure 1:** Contingency tables of KRAS mutations  
  **Script:** `data_analysis/enrichmentWeakStrongDrivers.R`  
  **Input data:**
  - `reduced_data/combinedKrasListDriverAnnotated.csv`
  - `reduced_data/cbioKrasDataNoPoleEdm.csv`
  - `reduced_data/cbioKrasDataWPoleEdm.csv`

- **Supplementary Figure 2:** Posterior of non-hypermutant weak driver rate  
  **Script:** `data_analysis/fig6DtoHPlotting/postDensMuWeak.R`  
  **Input data:** `reduced_data/simData/20260108.simKrasInf50Ksims.rds`

- **Supplementary Table 1:** Summary of studies used  
  **Script:** `data_analysis/annotatePublicMutData.R`  
  **Data:** `reduced_data/studySummaryTable.csv`

- **Supplementary Table 2:** Annotated KRAS driver classification  
  **Script:** `data_analysis/annotateKrasDriverClass.R`  
  **Data:** `reduced_data/combinedKrasListDriverAnnotated.csv`

---

