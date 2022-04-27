# Deciphering the regulatory code of histone modifications in plants
## Overview
Histone modifications play a critical role in the dynamic regulation of spatial-temporal gene expression. 
However, in plants, what genetic mechanisms determine the modification level of histone marks is still elusive. 
Using machine learning approach, we performed a systematic and comprehensive study on the cis-regulatory grammar of 19 histone marks in Arabidopsis thaliana, rice and maize.

## The schematic flowchart of GSminer <img src="Rcode/figures/flowchart.png" />





## Usage

### Predict histone modification levels
``` r
# The easiest way to predict histone modifications levels with DNA sequences using our R script:
source("R1_AthLeaf_predict_HMs.R")
peakKmerResReads.df <- getSeqFeature(peakFile=markers[i], k=6)
peakFeatures.df <- getMarkerReadsSignal(peakKmerResReads.df=peakKmerResReads.df,
                                             bamFile=bamFiles[i],peakFile = markers[i])
AthLeafPerformance <- cv(peakKmerRes.df = peakFeatures.df)
```

### Extracting mark motifs

``` r
require(gtools)
require(Biostrings)
source("R5_PWM_1201.R")
labc <-getActiveVar(peakSignal.df=mark_pos_neg_KmerRes.df)
AthLeafAllLASSOImp[[i]] <- labc
aa1 <- getEnrichmentImp(LASSO.impAll=labc,peakFile=markers[i])
AthLeafAllEnrichmentImp[[i]] <- aa1
bb1 <- getMismatch_k_mers_list(KmerWeight=aa1$newKmerWeight,enrichment.imp=aa1$EF1_df)
AthLeafAllMismatch_k_mers_list[[i]] <- bb1
markPWMs <- getEachMarkPWMs(kmerList = bb1,KmerWeight = aa1$newKmerWeight)
write_meme(markPWMs, file=markmotif.meme)
```

## Citation
Zhaohong Li#, Dongwei Li#, Ye Li, Xiaoping Guo, Ruolin Yang*: Deciphering the regulatory code of histone modifications in plants. 2022