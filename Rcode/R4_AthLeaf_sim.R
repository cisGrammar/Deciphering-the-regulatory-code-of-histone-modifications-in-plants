############################################################################################
# Date: 2019 10 16
# leaf in Arabidopsis thaliana
# simulation : shuffle the sequence; add the noise to the signal; sample the signal
# lidongwei@120.95.131.135 -p 2222
# scp yanglab@120.95.131.236:/home/yanglab/Allwork/part1/Rcommands/Ath/leaf/simulation/AthLeaf_simR.R /data1/lidongwei/ldw/data/Rcommands/Ath/leaf
############################################################################################

library(ranger)
library(rtracklayer)
library(Biostrings)
library(universalmotif)
library(gtools)
library(Rsamtools)
library(bamsignals)

setwd("/data1/lidongwei/ldw/data/Markers/AthNarrowpeak/leaf/")
markers <- list.files(pattern = "bam_peaks.narrowPeak")
bamFiles <- list.files(path = "/data1/lidongwei/ldw/data/Markers/AthNarrowpeak/leaf/bam/",
                       pattern = "bam$",full.names = T)


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


getSeqFeature <- function(peakSeq, k = 6)
{
  require(Biostrings)
  require(rtracklayer)
  freqMat <- oligonucleotideFrequency(peakSeq, width = k, step = 1)
  seqFeature <- collapseSeq(k = k)
  seqFeatureMat <-
    matrix(0,
           nrow = length(peakSeq),
           ncol = length(seqFeature))
  colnames(seqFeatureMat) <- seqFeature
  rownames(seqFeatureMat) <- names(peakSeq)
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

#### simulation_1 : shuffle the sequence
simulate_k <- function(peakFile,bamFile, keepK=2, k=6, niter=1000, AthFaFile="/data1/lidongwei/ldw/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
{
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakSeq <- getSeq(FaFile(AthFaFile), trim(peak))
  names(peakSeq) <- peak$name
  signal <-log2(bamCount(bamFile,peak,verbose=FALSE)+1)
  truePeakKmerRes.df <- getSeqFeature(peakSeq = peakSeq, k = k)
  truePeakKmerRes.df$signal <- signal
  trueModel <- ranger(signal~., truePeakKmerRes.df, num.threads = 40)
  trueRsquared <- trueModel$r.squared 
  truePred <- predict(trueModel, truePeakKmerRes.df)
  trueCor <- cor.test(truePred$prediction, truePeakKmerRes.df$signal)$estimate
  simRes <- data.frame(Rsquared=vector("numeric", length = niter), cor=vector("numeric", length = niter),
                       true_Rsquared= trueRsquared,true_cor =trueCor)
  for(i in 1:niter)
  {
    message(sprintf("in %s %d ...", peakFile, i))
    set.seed(i)
    peakSeqSim <- shuffle_sequences(peakSeq, k = keepK, method = "markov")
    peakKmerRes.df <- getSeqFeature(peakSeq = peakSeqSim, k = k)  
    peakKmerRes.df$signal <- signal
    model <- ranger(signal~., peakKmerRes.df, num.threads = 40)
    simRes[i,1] <- model$r.squared
    pred <- predict(model, peakKmerRes.df)
    cor <- cor.test(pred$prediction, peakKmerRes.df$signal)$estimate
    simRes[i,2] <- cor
    message(sprintf("in %s %d: sim_Rsquared: %f sim_cor: %f true_Rsquared: %f true_cor: %f ", peakFile, i, simRes[i,1], simRes[i,2], trueRsquared, trueCor))
  }
  save(simRes, file=paste0(unlist(strsplit(peakFile, split="\\."))[1], "_keep2_sim.RData"))
  simRes
}


#### simulation_2 : add the noise
noiseSignal <- function(peakFile,bamFile,k=6, niter=1000, AthFaFile="/data1/lidongwei/ldw/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
{
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakSeq <- getSeq(FaFile(AthFaFile), trim(peak))
  names(peakSeq) <- peak$name
  peakKmerRes.df <- getSeqFeature(peakSeq = peakSeq, k = k)
  signal <- log2(bamCount(bamFile,peak,verbose=FALSE)+1)
  peakKmerRes.df$signal <- signal
  trueModel <- ranger(signal~., peakKmerRes.df, num.threads = 40)
  trueRsquared <- trueModel$r.squared 
  truePred <- predict(trueModel, peakKmerRes.df)
  trueCor <- cor.test(truePred$prediction, peakKmerRes.df$signal)$estimate
  simRes <- data.frame(Rsquared=vector("numeric", length = niter), cor=vector("numeric", length = niter),
                       true_Rsquared= trueRsquared,true_cor =trueCor)
  for(i in 1:niter)
  {
    set.seed(i)
    noise <- rnorm(length(signal), mean = mean(signal), sd = sd(signal))
    peakKmerRes.df$signal <- signal + noise 
    model <- ranger(signal~., peakKmerRes.df, num.threads = 40)
    simRes[i,1] <- model$r.squared
    pred <- predict(model, peakKmerRes.df)
    cor <- cor.test(pred$prediction, signal)$estimate
    simRes[i,2] <- cor
    message(sprintf("in %s %d: sim_Rsquared: %f sim_cor: %f true_Rsquared: %f true_cor: %f ", peakFile, i, simRes[i,1], simRes[i,2], trueRsquared, trueCor))
  }
  save(simRes, file=paste0(unlist(strsplit(peakFile, split="\\."))[1], "_noiseSignal_sim.RData"))
  simRes
}

#### simulation_3 : sample the signal
sampleSignal <- function(peakFile,bamFile,k=6, niter=1000, AthFaFile="/data1/lidongwei/ldw/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
{
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakSeq <- getSeq(FaFile(AthFaFile), trim(peak))
  names(peakSeq) <- peak$name
  peakKmerRes.df <- getSeqFeature(peakSeq = peakSeq, k = k)
  signal <- log2(bamCount(bamFile,peak,verbose=FALSE)+1)
  peakKmerRes.df$signal <- signal
  trueModel <- ranger(signal~., peakKmerRes.df, num.threads = 40)
  trueRsquared <- trueModel$r.squared 
  truePred <- predict(trueModel, peakKmerRes.df)
  trueCor <- cor.test(truePred$prediction, peakKmerRes.df$signal)$estimate
  simRes <- data.frame(Rsquared=vector("numeric", length = niter), cor=vector("numeric", length = niter),
                       true_Rsquared= trueRsquared,true_cor =trueCor)
  for(i in 1:niter)
  {
    set.seed(i)
    peakKmerRes.df$signal <- sample(signal ,length(signal),replace = F) 
    model <- ranger(signal~., peakKmerRes.df, num.threads = 40)
    simRes[i,1] <- model$r.squared
    pred <- predict(model, peakKmerRes.df)
    cor <- cor.test(pred$prediction, signal)$estimate
    simRes[i,2] <- cor
    message(sprintf("in %s %d: sim_Rsquared: %f sim_cor: %f true_Rsquared: %f true_cor: %f ", peakFile, i, simRes[i,1], simRes[i,2], trueRsquared, trueCor))
  }
  save(simRes, file=paste0(unlist(strsplit(peakFile, split="\\."))[1], "_sampleSignal_sim.RData"))
  simRes
}

## results_1:
simMarkers_shuffle <-lapply(markers,function(x) simulate_k(peakFile = x,bamFile =bamFiles[which(markers==x)] ,niter = 1000))  
names(simMarkers_shuffle) <- gsub("\\.bam_peaks\\.narrowPeak","",markers)
save(simMarkers_shuffle, file="/data1/lidongwei/ldw/data/results/Ath/leaf/simulation/all_sim_keep2Signal.RData")

## result_2:
simMarkers_noise <-lapply(markers,function(x) noiseSignal(peakFile = x,bamFile =bamFiles[which(markers==x)] ,niter = 1000))  
names(simMarkers_noise) <- gsub("\\.bam_peaks\\.narrowPeak","",markers)
save(simMarkers_noise, file="/data1/lidongwei/ldw/data/results/Ath/leaf/simulation/all_sim_noiseSignal.RData")

## result_3:
simMarkers_sample <-lapply(markers,function(x) sampleSignal(peakFile = x,bamFile =bamFiles[which(markers==x)] ,niter = 1000))  
names(simMarkers_sample) <- gsub("\\.bam_peaks\\.narrowPeak","",markers)
save(simMarkers_sample, file="/data1/lidongwei/ldw/data/results/Ath/leaf/simulation/all_sim_sampleSignal.RData")

