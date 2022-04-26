############################################################################
# leaf in Arabidopsis thaliana 
# using the 6-mer under the peak to predict HMs
# model : linear regresison
#############################################################################
# leaf in Arabidopsis thaliana
setwd("/home/yanglab/data/Markers/AthNarrowpeak/leaf/new")   

markers <- list.files(pattern = "bam_peaks.narrowPeak")
bamFiles <- list.files(path = "/home/yanglab/data/Markers/AthNarrowpeak/leaf/new/bam/",
                       pattern = "bam$",full.names = T)

## get the sequence feature

collapseSeq <- function(k = 6,
                        nucleotides = c("A", "T", "C", "G"))
{
  require(gtools)
  require(Biostrings)
  allperms <-
    apply(permutations(length(nucleotides), k, nucleotides, repeats.allowed = T),
          1,
          paste0,
          collapse = "")
  res <- vector(mode = "character")
  for (i in seq_along(allperms))
  {
    revComp <- as.character(reverseComplement(DNAString(allperms[i])))
    if (revComp == allperms[i]) {
      res <- c(res, allperms[i])
      next
    } else{
      if (!revComp %in% res && !allperms[i] %in% res)
        res <- c(res, allperms[i])
    }
  }
  return(res)
}


getSeqFeature <- function(AthFaFile="/home/yanglab/lic/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                          peakFile, k = 6)
{
  require(Biostrings)
  require(rtracklayer)
  require(Rsamtools)
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakSeq <- getSeq(FaFile(AthFaFile), peak)
  freqMat <- oligonucleotideFrequency(peakSeq, width = k, step = 1)
  seqFeature <- collapseSeq(k = k)
  seqFeatureMat <-
    matrix(0,
           nrow = length(peakSeq),
           ncol = length(seqFeature))
  colnames(seqFeatureMat) <- seqFeature
  rownames(seqFeatureMat) <- peak$name
  for (kmer in colnames(seqFeatureMat))
  {
    revComp <- as.character(reverseComplement(DNAString(kmer)))
    if (kmer == revComp)
    {
      if (kmer %in% colnames(freqMat))
      {
        seqFeatureMat[, kmer] <- freqMat[, kmer]
      } else{
        seqFeatureMat[, kmer] <- 0
      }
      next
    }
    if (kmer != revComp)
    {
      if (kmer %in% colnames(freqMat) && revComp %in% colnames(freqMat))
      {
        seqFeatureMat[, kmer] <- freqMat[, kmer] + freqMat[, revComp]
      } else if (kmer %in% colnames(freqMat) &&
                 !revComp %in% colnames(freqMat))
      {
        seqFeatureMat[, kmer] <- freqMat[, kmer]
      } else if (!kmer %in% colnames(freqMat) &&
                 revComp %in% colnames(freqMat))
      {
        seqFeatureMat[, kmer] <- freqMat[, revComp]
      } else{
        seqFeatureMat[, kmer] <- 0
      }
      next
    }
  }
  peakKmerRes.df <- as.data.frame(seqFeatureMat)
  peakKmerRes.df
}

#### signal = reads
getMarkerReadsSignal <- function(peakKmerResReads.df, bamFile,peakFile)
{
  library(bamsignals)
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakKmerResReads.df$signal <- bamCount(bamFile,peak,verbose=FALSE)
  peakKmerResReads.df$signal <- log2(peakKmerResReads.df$signal +1)
  return(peakKmerResReads.df)
}

### to get the observed and predicted signal
cv_M_P <- function(peakKmerResReads.df, k=10)
{
  require(caret)
  folds <- createFolds(peakKmerResReads.df$signal, k = 10)
  predicted_ten_fold<-c()
  measured_ten_fold<-c()
  for (k in 1:length(folds))
  {
    cv_train <- peakKmerResReads.df[-folds[[k]],  ]
    cv_test <- peakKmerResReads.df[folds[[k]], ]
    cv_model <- lm(signal~., cv_train)
    cv_pred <- predict(cv_model, cv_test)
    predicted_ten_fold<-c(predicted_ten_fold,cv_pred)
    measured_ten_fold<-c(measured_ten_fold,cv_test$signal)
  }
  PredictPerformance <- data.frame(predicted = predicted_ten_fold,
                                   observed = measured_ten_fold)
  return(PredictPerformance)
}



AthLeafPerformances_obs_pre_lm <- vector("list", length = length(markers))
names(AthLeafPerformances_obs_pre_lm) <- markers



## get the performance of each mark

for (i in 1:length(bamFiles))
{
  peakKmerResReads.df <- getSeqFeature(peakFile=markers[i], k=6)
  peakKmerResReads.df <- getMarkerReadsSignal(peakKmerResReads.df=peakKmerResReads.df, 
                                              bamFile=bamFiles[i],peakFile = markers[i])
  AthLeafPerformances_obs_pre_lm[[i]] <- cv_M_P(peakKmerResReads.df = peakKmerResReads.df)
}

## save the data

saveRDS(AthLeafPerformances_obs_pre_lm,
        file = "/home/yanglab/Allwork/part1/results/Ath/leaf/lm/AthLeafPerformances_obs_pre_lm.RDS")


## plot the Scatter plot to show the observed signal and predicted signal

## exclude h2a 2021 7 10
# AthLeafPerformances_obs_pre_lm <- readRDS(file = "/home/yanglab/Allwork/part1/results/Ath/leaf/lm/AthLeafPerformances_obs_pre_lm.RDS")

AthLeafPerformances_obs_pre_lm <- AthLeafPerformances_obs_pre_lm[-1]
# list to data.frame

AthLeafPerformances_obs_pre_lm1 <- lapply(1:length(AthLeafPerformances_obs_pre_lm), function(i) {
  marker <- gsub("(\\.bam)?_peaks\\.narrowPeak","",names(AthLeafPerformances_obs_pre_lm)[i])
  AthLeafPerformances_obs_pre_lm[[i]]$Marker <- marker
  AthLeafPerformances_obs_pre_lm[[i]]
})

res <- data.frame()
for(i in 1:length(AthLeafPerformances_obs_pre_lm1))
  res <- rbind(res, AthLeafPerformances_obs_pre_lm1[[i]])


# plot 
p <- ggplot(data = res, aes(x = observed,y = predicted)) + 
  geom_point() +
  stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
  theme_bw() +
  theme(panel.grid=element_blank())  +
  geom_abline(intercept=0,slope=1) +
  scale_fill_viridis_c() + facet_wrap(.~Marker, nrow = 3, scales = "free") 

library(plyr)
df.cor <- ddply(res, .(Marker), function(val) sprintf("r==%.3f", cor(val$predicted, val$observed)))

pp <- p+geom_text(data=df.cor, aes(x=0, y=17, label=V1), parse=TRUE,hjust = 0,vjust=-0.8,size=4)+
  scale_y_continuous(limits=c(0,20))+
  scale_x_continuous(limits = c(0,20))

# save : /home/yanglab/Allwork/part1/results/Ath/leaf/lm/AthLeafPerformances_obs_pre_lm_sub.pdf
# 9 * 7





















