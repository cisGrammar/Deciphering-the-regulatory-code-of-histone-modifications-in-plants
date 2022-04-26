##################################################################################
# leaf in Arabidopsis thaliana
# divide the peak into bin(upstream & downstream 1000bp ) to predict
# Date : 2019 10 22
####################################################################################

setwd("/home/yanglab/data/Markers/AthNarrowpeak/leaf/new/")   

marks <- list.files(pattern = "bam_peaks.narrowPeak")
bamFiles <- list.files(path = "/home/yanglab/data/Markers/AthNarrowpeak/leaf/new/bam/",
                       pattern = "bam$",full.names = T)
summitFiles <- list.files(path = "/home/yanglab/data/Markers/AthNarrowpeak/leaf/new/peakSummit",
                          pattern = "summits.bed",full.names = T)


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

getBinSeqFeature <- function(peak.gr, k=6, AthFaFile)
{
  require(Biostrings)
  require(rtracklayer)
  require(Rsamtools)
  #peak <- import(peak.gr)
  peak <- peak.gr[seqnames(peak.gr) != "Mt" & seqnames(peak.gr) != "Pt"]
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

extend <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}


seqFeatureBinTile <- function(peakSummitFile, k=6, refFile, binLen = 100, expand = 10)
{
  require(Rsamtools)
  require(rtracklayer)
  require(GenomicRanges)
  ref.gr <- scanFaIndex(refFile, as=c("GRangesList", "GRanges"))
  peakSummit <- import(peakSummitFile)
  peakSummit <- peakSummit[seqnames(peakSummit) != "Mt" & seqnames(peakSummit) != "Pt"]
  seqlevels(peakSummit) <- seqlevels(ref.gr)
  seqlengths(peakSummit) <- seqlengths(ref.gr)
  peak_extend <- suppressWarnings(
    trim(extend(peakSummit, upstream = expand * binLen, downstream = expand * binLen))
  )
  peakTile.list <- tile(peak_extend, n=2*expand)
  
  names(peakTile.list) <- peakSummit$name
  return(peakTile.list)
}

cv <- function(PeakKmerRes.df, k=10)
{
  require(caret)
  folds <- createFolds(PeakKmerRes.df$signal, k = k)
  lapply(folds, function(x)
  {
    cv_train <- PeakKmerRes.df[-x, ]
    cv_test <- PeakKmerRes.df[x, ]
    cv_model <- ranger(signal~., cv_train,num.threads = 7)
    cv_pred <- predict(cv_model, cv_test)
    corRes <- cor.test(cv_pred$prediction,cv_test$signal)
    res <- c(corRes$p.value, corRes$estimate, cv_model$r.squared)
    names(res) <- c("pvalue", "cor", "rsquared")
    return(res)
  }
  )
}


seqBinTraining <- function(peakFile, peakTile.list, refFile,bamFile, k=6)
{
  require(ranger)
  require(S4Vectors)
  library(bamsignals)
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  if(length(peak) != length(peakTile.list))
    stop("peak file length must equal with the length of tiled list!")
  lapply(peak$name, function(id)
  {
    message(id)
    tmp <- peak[peak$name == id]
    peakTile.list[[id]]$signal <<- log2(bamCount(bamFile,tmp,verbose=FALSE) +1)
  })
  
  res <- vector("list", length = length(peakTile.list[[1]]))
  for(i in 1:length(peakTile.list[[1]]))
  {
    message(format(Sys.time(), "%a %b %d %X %Y")," 10 fold cv in bin: ",i)
    gr <- unlist(as(lapply(peakTile.list, function(x) `[`(x,i)), "GRangesList"))
    binKmer.df <- getBinSeqFeature(peak.gr = gr, k = k, AthFaFile = refFile)
    binKmer.df$signal <- gr$signal
    res[[i]] <- as.data.frame(cv(binKmer.df))
  }
  res
}


plotBin <- function(performance)
{
  performance <- lapply(performance, function(x){ 
    nfold <- colnames(x)
    res <- as.data.frame(t(x))
    res$nfold <- nfold
    res
  })
  lapply(1:length(performance), function(x) {
    performance[[x]]$bin <<- x
  })
  df <- do.call("rbind", performance)
  require(ggplot2)
  ggplot(data = df, aes(x=bin, y=cor, group=bin)) + geom_boxplot()
}

AthLeafPerformance_bin <- vector("list", length = length(marks))
names(AthLeafPerformance_bin) <- marks



for (i in 1:length(marks))
{
  peakTile.list <- seqFeatureBinTile(peakSummitFile = summitFiles[i],
                                     refFile = "/home/yanglab/lic/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
  AthLeafPerformance_bin[[i]]<-seqBinTraining(peakFile=marks[i], 
                                                  peakTile.list = peakTile.list, 
                                                  refFile = "/home/yanglab/lic/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                                                  bamFile = bamFiles[i])
  saveRDS(AthLeafPerformance_bin[[i]],
          file=paste0(unlist(strsplit(marks[i], split="\\."))[1], "_bin.RData"))
}

saveRDS(AthLeafPerformance_bin,"AthLeafPerformance_bin")


