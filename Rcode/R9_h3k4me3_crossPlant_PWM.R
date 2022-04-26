##############################################################################
# This part tries to explore weather sequences could predict histone marks across
#  plants, which to explore the degree of evolutionary conservation 
# the relationship between motifs and histone marks.
# Data : 2021 6 21
# pwd : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction
##############################################################################
#########

# 1 data 
 ## pwd :  /home/yanglab/Allwork/part4/data/Plant_H3K4me3
AthH3K4me3Peak <- "/home/yanglab/Allwork/part4/data/Plant_H3K4me3/AthLeafH3K4me3.bam_peaks.narrowPeak"
AthH3K4me3Bam <- "/home/yanglab/Allwork/part4/data/Plant_H3K4me3/AthLeafH3K4me3.bam"
AthH3K4me3PWM <- "/home/yanglab/Allwork/reWork/part3/cross_plant/compare_PWMs/AthLeafH3K4me3.meme"
# AthFaFile = "/home/yanglab/lic/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

RiceH3K4me3Peak <- "/home/yanglab/Allwork/part4/data/Plant_H3K4me3/RiceLeafH3K4me3.bam_peaks.narrowPeak"
RiceH3K4me3Bam <- "/home/yanglab/Allwork/part4/data/Plant_H3K4me3/RiceLeafH3K4me3.bam"
RiceH3K4me3PWM <- "/home/yanglab/Allwork/reWork/part3/cross_plant/compare_PWMs/RiceLeafH3K4me3.meme"
# RiceFaFile <- "/home/yanglab/lic/GenomeRef/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"

MaizeH3k4me3Peak <- "/home/yanglab/Allwork/part4/data/Plant_H3K4me3/MaizeLeafH3K4me3.bam_peaks.narrowPeak" 
MaizeH3k4me3Bam <- "/home/yanglab/Allwork/part4/data/Plant_H3K4me3/MaizeLeafH3K4me3.bam"
MaizeH3k4me3PWM <- "/home/yanglab/Allwork/reWork/part3/cross_plant/compare_PWMs/MaizeLeafH3K4me3.meme"
# MaizeFaFile <- "/home/yanglab/lic/GenomeRef/Zea_mays.AGPv4.dna.toplevel.fa"

############################
## Ath
getPWMPeakCount_Ath<- function(markPWM,peakFile)
{
  library(Rsamtools)
  library(Biostrings)
  library(rtracklayer)
  faFile="/home/yanglab/lic/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakSeq <- getSeq(FaFile(faFile), peak)
  library(universalmotif)
  markPWM <- read_meme(markPWM) 
  PeakmarkPWMCounts1 <- matrix(ncol = length(markPWM),
                               nrow = length(peakSeq))
  motifNames <-vector()
  for (i in 1:length(markPWM))
    motifNames[i] <- markPWM[[i]]@name
  for (i in 1:length(markPWM)){
    PeakmarkPWMCounts1[,i] <- as.numeric(lapply(peakSeq, function(s) countPWM(markPWM[[i]]@motif, s, min.score="85%") ))
  }
  
  PeakmarkPWMCounts2 <- matrix(ncol = length(markPWM),
                               nrow = length(peakSeq))
  for (i in 1:length(markPWM)){
    PeakmarkPWMCounts2[,i] <- as.numeric(lapply(peakSeq, function(s) countPWM(reverseComplement(markPWM[[i]]@motif), s, min.score="85%") ))
  }
  PeakmarkPWMCounts <-as.data.frame(PeakmarkPWMCounts1 + PeakmarkPWMCounts2)
  colnames(PeakmarkPWMCounts) <- motifNames
  rownames(PeakmarkPWMCounts) <- paste0("Peak","_",1:length(peak))
  return(PeakmarkPWMCounts)
}


#### signal = reads
getMarkerReadsSignal_Ath <- function(peakPWMResReads.df, bamFile,peakFile)
{
  library(bamsignals)
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakPWMResReads.df$signal <- bamCount(bamFile,peak,verbose=FALSE)
  peakPWMResReads.df$signal <- log2(peakPWMResReads.df$signal +1)
  return(peakPWMResReads.df)
}


## Rice

getPWMPeakCount_Rice<- function(markPWM,peakFile)
{
  library(Rsamtools)
  library(Biostrings)
  library(rtracklayer)
  faFile="/home/yanglab/lic/GenomeRef/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) %in% as.character(c(1:12))]
  peakSeq <- getSeq(FaFile(faFile), peak)
  library(universalmotif)
  markPWM <- read_meme(markPWM) 
  PeakmarkPWMCounts1 <- matrix(ncol = length(markPWM),
                               nrow = length(peakSeq))
  motifNames <-vector()
  for (i in 1:length(markPWM))
    motifNames[i] <- markPWM[[i]]@name
  for (i in 1:length(markPWM)){
    PeakmarkPWMCounts1[,i] <- as.numeric(lapply(peakSeq, function(s) countPWM(markPWM[[i]]@motif, s, min.score="85%") ))
  }
  
  PeakmarkPWMCounts2 <- matrix(ncol = length(markPWM),
                               nrow = length(peakSeq))
  for (i in 1:length(markPWM)){
    PeakmarkPWMCounts2[,i] <- as.numeric(lapply(peakSeq, function(s) countPWM(reverseComplement(markPWM[[i]]@motif), s, min.score="85%") ))
  }
  PeakmarkPWMCounts <-as.data.frame(PeakmarkPWMCounts1 + PeakmarkPWMCounts2)
  colnames(PeakmarkPWMCounts) <- motifNames
  rownames(PeakmarkPWMCounts) <- paste0("Peak","_",1:length(peak))
  return(PeakmarkPWMCounts)
}


#### signal = reads
getMarkerReadsSignal_Rice <- function(peakPWMResReads.df, bamFile,peakFile)
{
  library(bamsignals)
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) %in% as.character(c(1:12))]
  peakPWMResReads.df$signal <- bamCount(bamFile,peak,verbose=FALSE)
  peakPWMResReads.df$signal <- log2(peakPWMResReads.df$signal +1)
  return(peakPWMResReads.df)
}


## Maize

getPWMPeakCount_Maize<- function(markPWM,peakFile)
{
  library(Rsamtools)
  library(Biostrings)
  library(rtracklayer)
  faFile="/home/yanglab/lic/GenomeRef/Zea_mays.AGPv4.dna.toplevel.fa"
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) %in% as.character(c(1:10))]
  peakSeq <- getSeq(FaFile(faFile), peak)
  library(universalmotif)
  markPWM <- read_meme(markPWM) 
  PeakmarkPWMCounts1 <- matrix(ncol = length(markPWM),
                               nrow = length(peakSeq))
  motifNames <-vector()
  for (i in 1:length(markPWM))
    motifNames[i] <- markPWM[[i]]@name
  for (i in 1:length(markPWM)){
    PeakmarkPWMCounts1[,i] <- as.numeric(lapply(peakSeq, function(s) countPWM(markPWM[[i]]@motif, s, min.score="85%") ))
  }
  
  PeakmarkPWMCounts2 <- matrix(ncol = length(markPWM),
                               nrow = length(peakSeq))
  for (i in 1:length(markPWM)){
    PeakmarkPWMCounts2[,i] <- as.numeric(lapply(peakSeq, function(s) countPWM(reverseComplement(markPWM[[i]]@motif), s, min.score="85%") ))
  }
  PeakmarkPWMCounts <-as.data.frame(PeakmarkPWMCounts1 + PeakmarkPWMCounts2)
  colnames(PeakmarkPWMCounts) <- motifNames
  rownames(PeakmarkPWMCounts) <- paste0("Peak","_",1:length(peak))
  return(PeakmarkPWMCounts)
}


#### signal = reads
getMarkerReadsSignal_Maize <- function(peakPWMResReads.df, bamFile,peakFile)
{
  library(bamsignals)
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) %in% as.character(c(1:10))]
  peakPWMResReads.df$signal <- bamCount(bamFile,peak,verbose=FALSE)
  peakPWMResReads.df$signal <- log2(peakPWMResReads.df$signal +1)
  return(peakPWMResReads.df)
}


## model prediction

cv <- function(peakPWMRes.df, k=10)
{
  require(ranger)
  require(caret)
  folds <- createFolds(peakPWMRes.df$signal, k = k)
  lapply(folds, function(x)
  {
    cv_train <- peakPWMRes.df[-x, ]
    cv_test <- peakPWMRes.df[x, ]
    cv_model <- ranger(signal~., cv_train,num.threads = 6)
    cv_pred <- predict(cv_model, cv_test)
    corRes <- cor.test(cv_pred$prediction,cv_test$signal)
    return(c(corRes$p.value, corRes$estimate, cv_model$r.squared))
  }
  )
}


acrossCV <- function(peakKmerResTrain.df,peakKmerResTest.df, niter=10)
{
  require(ranger)
  require(caret)
  cv_train<- peakKmerResTrain.df
  cv_model <- ranger(signal~., cv_train,num.threads = 6)
  simRes <- data.frame(p_values=vector("numeric", length = niter), cor=vector("numeric", length = niter))
  for(i in 1:niter)
  {
    set.seed(i)
    indexs <- sample(1:nrow(peakKmerResTest.df),nrow(peakKmerResTest.df)*0.8)
    peakKmerResTest_new.df <- peakKmerResTest.df[indexs,]
    PredPromance <- predict(cv_model,peakKmerResTest_new.df)
    corRes <- cor.test(PredPromance$prediction,peakKmerResTest_new.df$signal)
    simRes[i,1] <- corRes$p.value
    simRes[i,2] <- corRes$estimate
  }
  return(simRes)
}



# 2 prediction

############### H3K4me3 ################################################################

AthLeafH3K4me3PeakKmerRes_AthPwWM.df <- getPWMPeakCount_Ath(markPWM = AthH3K4me3PWM,
                                                peakFile = AthH3K4me3Peak)
AthLeafH3K4me3PeakKmerRes_AthPwWM.df <- getMarkerReadsSignal_Ath(AthLeafH3K4me3PeakKmerRes_AthPwWM.df, 
                                                        bamFile = AthH3K4me3Bam,
                                                        peakFile=AthH3K4me3Peak)
saveRDS(AthLeafH3K4me3PeakKmerRes_AthPwWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/AthH3K4me3PWMCount_df_AthPwWM.RDS")


RiceLeafH3K4me3PeakKmerRes_AthPwWM.df <- getPWMPeakCount_Rice(markPWM = AthH3K4me3PWM,
                                                     peakFile = RiceH3K4me3Peak)
RiceLeafH3K4me3PeakKmerRes_AthPwWM.df <- getMarkerReadsSignal_Rice(RiceLeafH3K4me3PeakKmerRes_AthPwWM.df, 
                                                          bamFile = RiceH3K4me3Bam,
                                                          peakFile=RiceH3K4me3Peak)
saveRDS(RiceLeafH3K4me3PeakKmerRes_AthPwWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/RiceH3K4me3PWMCount_df_AthPwWM.RDS")


MaizeLeafH3K4me3PeakKmerRes_AthPwWM.df <- getPWMPeakCount_Maize(markPWM = AthH3K4me3PWM,
                                                      peakFile = MaizeH3k4me3Peak)
MaizeLeafH3K4me3PeakKmerRes_AthPwWM.df <- getMarkerReadsSignal_Maize(MaizeLeafH3K4me3PeakKmerRes_AthPwWM.df, 
                                                           bamFile = MaizeH3k4me3Bam,
                                                           peakFile=MaizeH3k4me3Peak)
saveRDS(MaizeLeafH3K4me3PeakKmerRes_AthPwWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/MaizeH3K4me3PWMCount_df_AthPwWM.RDS")




AthLeafH3K4me3CV_result_AthPwWM <- cv(AthLeafH3K4me3PeakKmerRes_AthPwWM.df)
AthLeafH3K4me3CV_result_AthPwWM <- as.data.frame(do.call("rbind", AthLeafH3K4me3CV_result_AthPwWM))


## Ath to Rice & Maize

AthLeaf_RiceLeafH3K4me3CV_result_AthPwWM <- acrossCV(AthLeafH3K4me3PeakKmerRes_AthPwWM.df,RiceLeafH3K4me3PeakKmerRes_AthPwWM.df)
#AthLeaf_RiceLeafH3K4me3CV_result <- as.data.frame(do.call("rbind", AthLeaf_RiceLeafH3K4me3CV_result))

AthLeaf_MaizeLeafH3K4me3CV_result_AthPwWM <- acrossCV(AthLeafH3K4me3PeakKmerRes_AthPwWM.df,MaizeLeafH3K4me3PeakKmerRes_AthPwWM.df)
#AthLeaf_MaizeLeafH3K4me3CV_result <- as.data.frame(do.call("rbind", AthLeaf_MaizeLeafH3K4me3CV_result))


Ath_Rice_Maize_AthPwWM <- vector("list")
Ath_Rice_Maize_AthPwWM$AthLeafH3K4me3CV_result_AthPwWM <-AthLeafH3K4me3CV_result_AthPwWM
Ath_Rice_Maize_AthPwWM$AthLeaf_RiceLeafH3K4me3CV_result_AthPwWM <-AthLeaf_RiceLeafH3K4me3CV_result_AthPwWM
Ath_Rice_Maize_AthPwWM$AthLeaf_MaizeLeafH3K4me3CV_result_AthPwWM <- AthLeaf_MaizeLeafH3K4me3CV_result_AthPwWM

saveRDS(Ath_Rice_Maize_AthPwWM,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df//Ath_Rice_Maize_H3K4me3_results")




aa1 <- lapply(Ath_Rice_Maize_AthPwWM,function(x) mean(x[,2]))




############### H3K4me3 PWM of Rice ################################################################

AthLeafH3K4me3PeakKmerRes_RicePwWM.df <- getPWMPeakCount_Ath(markPWM = RiceH3K4me3PWM,
                                                    peakFile = AthH3K4me3Peak)
AthLeafH3K4me3PeakKmerRes_RicePwWM.df <- getMarkerReadsSignal_Ath(AthLeafH3K4me3PeakKmerRes_RicePwWM.df, 
                                                         bamFile = AthH3K4me3Bam,
                                                         peakFile=AthH3K4me3Peak)
saveRDS(AthLeafH3K4me3PeakKmerRes_RicePwWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/AthH3K4me3PWMCount_df_RicePwWM.RDS")


RiceLeafH3K4me3PeakKmerRes_RicePwWM.df <- getPWMPeakCount_Rice(markPWM = RiceH3K4me3PWM,
                                                      peakFile = RiceH3K4me3Peak)
RiceLeafH3K4me3PeakKmerRes_RicePwWM.df <- getMarkerReadsSignal_Rice(RiceLeafH3K4me3PeakKmerRes_RicePwWM.df, 
                                                           bamFile = RiceH3K4me3Bam,
                                                           peakFile=RiceH3K4me3Peak)
saveRDS(RiceLeafH3K4me3PeakKmerRes.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/RiceH3K4me3PWMCount_df_RicePwWM.RDS")


MaizeLeafH3K4me3PeakKmerRes_RicePwWM.df <- getPWMPeakCount_Maize(markPWM = RiceH3K4me3PWM,
                                                        peakFile = MaizeH3k4me3Peak)
MaizeLeafH3K4me3PeakKmerRes_RicePwWM.df <- getMarkerReadsSignal_Maize(MaizeLeafH3K4me3PeakKmerRes_RicePwWM.df, 
                                                             bamFile = MaizeH3k4me3Bam,
                                                             peakFile=MaizeH3k4me3Peak)
saveRDS(MaizeLeafH3K4me3PeakKmerRes_RicePwWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/MaizeH3K4me3PWMCount_df_RicePwWM.RDS")





RiceLeafH3K4me3CV_result_RicePwWM <- cv(RiceLeafH3K4me3PeakKmerRes_RicePwWM.df)
RiceLeafH3K4me3CV_result_RicePwWM <- as.data.frame(do.call("rbind", RiceLeafH3K4me3CV_result_RicePwWM))



RiceLeaf_AthLeafH3K4me3CV_result_RicePwWM <- acrossCV(RiceLeafH3K4me3PeakKmerRes.df,AthLeafH3K4me3PeakKmerRes_RicePwWM.df)
#RiceLeaf_AthLeafH3K4me3CV_result <- as.data.frame(do.call("rbind", RiceLeaf_AthLeafH3K4me3CV_result))

RiceLeaf_MaizeLeafH3K4me3CV_result_RicePwWM <- acrossCV(RiceLeafH3K4me3PeakKmerRes.df,MaizeLeafH3K4me3PeakKmerRes_RicePwWM.df)
#RiceLeaf_MaizeLeafH3K4me3CV_result <- as.data.frame(do.call("rbind", RiceLeaf_MaizeLeafH3K4me3CV_result))


Ath_Rice_Maize_RicePwWM <- vector("list")
Ath_Rice_Maize_RicePwWM$RiceLeafH3K4me3CV_result_RicePwWM <- RiceLeafH3K4me3CV_result_RicePwWM
Ath_Rice_Maize_RicePwWM$RiceLeaf_AthLeafH3K4me3CV_result_RicePwWM <- RiceLeaf_AthLeafH3K4me3CV_result_RicePwWM
Ath_Rice_Maize_RicePwWM$RiceLeaf_MaizeLeafH3K4me3CV_result_RicePwWM <- RiceLeaf_MaizeLeafH3K4me3CV_result_RicePwWM

saveRDS(Ath_Rice_Maize_RicePwWM,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/Ath_Rice_Maize_H3K4me3_results_RicePwWM.RDS")


aa2 <- lapply(Ath_Rice_Maize_RicePwWM,function(x) mean(x[,2]))



############### H3K4me3  PWM of Maize################################################################

AthLeafH3K4me3PeakKmerRes_MaizePWM.df <- getPWMPeakCount_Ath(markPWM = MaizeH3k4me3PWM,
                                                             peakFile = AthH3K4me3Peak)
AthLeafH3K4me3PeakKmerRes_MaizePWM.df <- getMarkerReadsSignal_Ath(AthLeafH3K4me3PeakKmerRes_MaizePWM.df, 
                                                                  bamFile = AthH3K4me3Bam,
                                                                  peakFile=AthH3K4me3Peak)
saveRDS(AthLeafH3K4me3PeakKmerRes_MaizePWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/AthH3K4me3PWMCount_df_MaizePWM.RDS")


RiceLeafH3K4me3PeakKmerRes_MaizePWM.df <- getPWMPeakCount_Rice(markPWM = MaizeH3k4me3PWM,
                                                               peakFile = RiceH3K4me3Peak)
RiceLeafH3K4me3PeakKmerRes_MaizePWM.df <- getMarkerReadsSignal_Rice(RiceLeafH3K4me3PeakKmerRes_MaizePWM.df, 
                                                                    bamFile = RiceH3K4me3Bam,
                                                                    peakFile=RiceH3K4me3Peak)
saveRDS(RiceLeafH3K4me3PeakKmerRes_MaizePWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/RiceH3K4me3PWMCount_df_MaizePWM.RDS")


MaizeLeafH3K4me3PeakKmerRes_MaizePWM.df <- getPWMPeakCount_Maize(markPWM = MaizeH3k4me3PWM,
                                                                 peakFile = MaizeH3k4me3Peak)
MaizeLeafH3K4me3PeakKmerRes_MaizePWM.df <- getMarkerReadsSignal_Maize(MaizeLeafH3K4me3PeakKmerRes_MaizePWM.df, 
                                                                      bamFile = MaizeH3k4me3Bam,
                                                                      peakFile=MaizeH3k4me3Peak)
saveRDS(MaizeLeafH3K4me3PeakKmerRes_MaizePWM.df,file = "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/MaizeH3K4me3PWMCount_df_MaizePWM.RDS")





MaizeLeafH3K4me3CV_result_MaizePWM <- cv(MaizeLeafH3K4me3PeakKmerRes_MaizePWM.df)
MaizeLeafH3K4me3CV_result_MaizePWM <- as.data.frame(do.call("rbind", MaizeLeafH3K4me3CV_result_MaizePWM))

## Maize to Ath & Rice

MaizeLeaf_AthLeafH3K4me3CV_result_MaizePWM <- acrossCV(MaizeLeafH3K4me3PeakKmerRes_MaizePWM.df,AthLeafH3K4me3PeakKmerRes_MaizePWM.df)
#MaizeLeaf_AthLeafH3K4me3CV_result <- as.data.frame(do.call("rbind", MaizeLeaf_AthLeafH3K4me3CV_result))

MaizeLeaf_RiceLeafH3K4me3CV_result_MaizePWM <- acrossCV(MaizeLeafH3K4me3PeakKmerRes_MaizePWM.df,RiceLeafH3K4me3PeakKmerRes_MaizePWM.df)
#MaizeLeaf_RiceLeafH3K4me3CV_result <- as.data.frame(do.call("rbind", MaizeLeaf_RiceLeafH3K4me3CV_result))




Ath_Rice_Maize_MaizePWM <- vector("list")
Ath_Rice_Maize_MaizePWM$MaizeLeafH3K4me3CV_result_MaizePWM <- MaizeLeafH3K4me3CV_result_MaizePWM
Ath_Rice_Maize_MaizePWM$MaizeLeaf_AthLeafH3K4me3CV_result_MaizePWM <- MaizeLeaf_AthLeafH3K4me3CV_result_MaizePWM
Ath_Rice_Maize_MaizePWM$MaizeLeaf_RiceLeafH3K4me3CV_result_MaizePWM <- MaizeLeaf_RiceLeafH3K4me3CV_result_MaizePWM

saveRDS(Ath_Rice_Maize_MaizePWM,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/crossPlant_PWMprediction/AthRiceMaizeH3K4me3PWMCount_df/Ath_Rice_Maize_H3K4me3_results_MaizePWM.RDS")




aa3 <- lapply(Ath_Rice_Maize_MaizePWM,function(x) mean(x[,2]))


aa3





