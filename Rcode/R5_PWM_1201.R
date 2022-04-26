########################################################################################
setwd("/home/yanglab/Allwork/reWork/part2/AthLeafPWMs/postivePeaks")   
markers <- list.files(pattern = "bam_peaks.narrowPeak")
markersNeg <- list.files(path = "/home/yanglab/data/Markers/AthNarrowpeak/leaf/new/negativePeak/",
                         pattern = "*.bed",full.names = T)
bamFiles <- list.files(path = "/home/yanglab/data/Markers/AthNarrowpeak/leaf/new/bam/",
                       pattern = "bam$",full.names = T)
# step_1 :
# get the dataset that contain the peaks and the negative dataset without any modifications of the specific mark

getDateSet <- function(peak.postive,peak.negtive){
  # Get the dataset contain the peak and negative peak.  
  #  
  # Args:  
  #   peak.postive : The peak file of the mark.  
  #   peak.negtive : The peak file of the region that without any modification of the specific mark.  
  # 
  # Returns:  
  #   The dataset contain postive and negative dataset.
  
  # load the packages
  
  require(rtracklayer)
  require(BiocGenerics)
  
  peak_postive <- import(peak.postive)
  peak_negtive <- import(peak.negtive)
  peak_postive <- peak_postive[seqnames(peak_postive) != "Mt" & seqnames(peak_postive) != "Pt"]
  peak_negtive <- peak_negtive[seqnames(peak_negtive) != "Mt" & seqnames(peak_negtive) != "Pt"]
  peak_negtive$name <- paste0("peak_neg","_",1:length(peak_negtive))
  peak_postive_gr <- granges(peak_postive)
  peak_postive_gr$name <- peak_postive$name
  peak_Data_set <- c(peak_postive_gr,peak_negtive)
  peak_Data_set <- BiocGenerics::sort(peak_Data_set)
  return(peak_Data_set)
}

 # mark_pos_neg_KmerRes.df <- getDateSet(markers[11],markersNeg[11])
# step_2 : get the k-mers data.frame

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
  require(Rsamtools)
  require(rtracklayer)
  peak <- peakFile
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

getMarkerReadsSignal <- function(peakKmerResReads.df, bamFile,peakFile)
{
  library(bamsignals)
  peak <-peakFile
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakKmerResReads.df$signal <- bamCount(bamFile,peak,verbose=FALSE)
  peakKmerResReads.df$signal <- log2(peakKmerResReads.df$signal +1)
  return(peakKmerResReads.df)
}


################################################################################

getActiveVar <- function(peakSignal.df)
{
  require(doMC)
  registerDoMC(cores=4)
  library(glmnet)
  x <- model.matrix(signal~.,peakSignal.df)[,-1]
  y <- peakSignal.df$signal
  cv.out <- cv.glmnet(x,y,nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = 1)
  bestlam <- cv.out$lambda.min
  co<-coef(cv.out, s = "lambda.min")
  co <- co[-1,]
  co1 <- as.matrix(co)
  colnames(co1) <- "weigth"
  Active.Index<-which(co1 > 0)
  variables<-names(co)[Active.Index]
  co1Varables <- subset(co1,co1[,1] >0)
  labc <- list(a=variables,b=co1Varables)
  return(labc)
}

########################################################################
###########################################################
################################################################################


#labc <-getActiveVar(peakSignal.df =mark_pos_neg_KmerRes.df)
#############################################################################


getEnrichmentImp<-function(LASSO.impAll, peakFile)
{
  library(rtracklayer)
  require(Rsamtools)
  library(Biostrings)
  require(universalmotif)
  library(dplyr)
  variables <-LASSO.impAll[[1]]
  co1Varables <- LASSO.impAll[[2]]
  faFile="/home/yanglab/lic/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
  peak <- import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakSeq <- getSeq(FaFile(faFile),peak)
  EF <-as.data.frame(matrix(numeric(0),ncol = 12,nrow = length(variables)))
  rownames(EF) <- variables
  colnames(EF) <- c("peakCounts1","peakCounts2","peakCounts",
                    "wgCounts1","wgCounts2","wgCounts",
                    "shCounts1","shCounts2","shCounts",
                    "EFwg","EFsh","CEF")
  ## peakCounts
  for (i in 1:length(variables))
    EF$peakCounts1[i] <- sum(vcountPattern(DNAString(variables[i]),peakSeq))
  ## peakCounts
  for (i in 1:length(variables))
    EF$peakCounts2[i] <- sum(vcountPattern(reverseComplement(DNAString(variables[i])),peakSeq))
  EF$peakCounts <- EF$peakCounts1 + EF$peakCounts2
  ### EFwg
  faFileSeq <-readDNAStringSet(faFile)
  faFileSeq <- faFileSeq[c(-6,-7),]
  for (i in 1:length(variables))
    EF$wgCounts1[i] <- sum(vcountPattern(DNAString(variables[i]),faFileSeq))
  for (i in 1:length(variables))
    EF$wgCounts2[i] <- sum(vcountPattern(reverseComplement(DNAString(variables[i])),faFileSeq))
  EF$wgCounts <- EF$wgCounts1 + EF$wgCounts2
  ## wshcounts
  peakSeqSh <- shuffle_sequences(peakSeq, k = 1, method = "euler")
  for (i in 1:length(variables))
    EF$shCounts1[i]<- sum(vcountPattern(DNAString(variables[i]),peakSeqSh))
  for (i in 1:length(variables))
    EF$shCounts2[i]<- sum(vcountPattern(reverseComplement(DNAString(variables[i])),peakSeqSh))
  EF$shCounts <- EF$shCounts1 + EF$shCounts2
  
  
  
  for (i in 1:length(variables))
    EF$EFsh[i] <- (EF$peakCounts[i]/sum(width(peak)))/(EF$shCounts[i]/sum(width(peak)))
  
  for (i in 1:length(variables))
    EF$EFwg[i] <- (EF$peakCounts[i]/sum(width(peak)))/(EF$wgCounts[i]/sum(width(faFileSeq)))
  
  
  for (i in 1:nrow(EF))
    EF$CEF[i] <- (EF$peakCounts[i]/sum(width(peak)))/(((EF$wgCounts[i]/sum(width(faFileSeq)))+
                                                         (EF$shCounts[i]/sum(width(peak))))*0.5)
  EF1 <- subset(EF,EF$EFwg >1 & EF$EFsh >1)
  CEF_weight <- EF %>% select(CEF)
  newKmerWeight <- co1Varables * CEF_weight
  colnames(newKmerWeight) <- "weight"
  newKmerWeight <- as.matrix(newKmerWeight)
  EFlist <- list(EF1_df = EF1,newKmerWeight = newKmerWeight)
  return(EFlist)
}

#aa1 <- getEnrichmentImp(LASSO.impAll = labc,peakFile = markers[11])


getMismatch_k_mers_list <- function(KmerWeight,enrichment.imp)
{
  aa <- data.frame(KmerWeight)
  aa$feature <- rownames(KmerWeight)
  bb<- data.frame(feature=rownames(enrichment.imp))
  bb$feature <- as.character(bb$feature)
  bb$weight <- aa$weight[(match(bb$feature,aa$feature))]
  bb1 <-bb[with(bb,order(weight,decreasing = T)),]
  bb2 <- bb1$feature
  
  strDis <- function(a, b)
  {
    dis <- 0    
    for(i in 1:nchar(a))
    {
      if(substr(a,i,i) != substr(b,i,i))
        dis <- dis + 1
    }
    return(dis)
  }
  
  strDisArr <- function(a, arr)
  {
    sapply(arr, function(x)strDis(a, x))
  }
  
  l=1
  allMismatch <- vector("list")
  while(length(bb2) >0)
  {
    bb11 <- names(strDisArr(bb2[1],aa$feature)[strDisArr(bb2[1],aa$feature) < 2])
    allMismatch[[l]] <- bb11
    bb2 <- setdiff(bb2,bb11)
    l = l+1
  }
  return(allMismatch)
}

#bb1 <- getMismatch_k_mers_list(KmerWeight = aa1$newKmerWeight,enrichment.imp = aa1$EF1_df)

pwm <- function(mismatch_k_mers,KmerWeight)
{
  library(universalmotif)
  ### pwm matrix
  pwmMatrix<- matrix(ncol=1,nrow=length(mismatch_k_mers))
  rownames(pwmMatrix) <- mismatch_k_mers
  for (i in rownames(pwmMatrix))
    pwmMatrix[i,1] <- KmerWeight[rownames(KmerWeight)==i]
  
  
  freqs <- matrix(ncol = 6,nrow = length(mismatch_k_mers))
  for (j in 1:6)
  {
    for (i in 1:length(mismatch_k_mers))
    {
      freqs[i,j] <- pwmMatrix[i]
    }
  }
  
  chrs <- matrix(ncol = 6,nrow = length(mismatch_k_mers))
  for (j in 1:6)
  {
    for (i in 1:length(mismatch_k_mers))
    {
      chrs[i,j] <- substring(mismatch_k_mers[i],1:6,1:6)[j]
    }
  }
  
  dna <-c("A","C","G","T")
  
  mats <- matrix(nrow = 4,ncol = 6)
  
  for (i in 1:length(dna))
  {
    for (j in 1:6)
    {
      indexs <- which(chrs[,j] == dna[i])
      mats[i,j] <- sum(freqs[,j][indexs])
    }
  }
  for (l1 in 1:nrow(mats)){
    for (l2 in 1:ncol(mats)){
      mats[l1,l2] <- mats[l1,l2] / sum(mat[,l2])
    }
  }
  kmers_pwm <- create_motif(mats,alphabet = "DNA", type = "PWM")
  return(kmers_pwm)
}

getEachMarkPWMs <- function(kmerList,KmerWeight)
{
  markPWMs <- vector("list",length = length(kmerList))
  for (l1 in 1:length(kmerList))
    markPWMs[[l1]] <- pwm(mismatch_k_mers = kmerList[[l1]],KmerWeight = aa1$newKmerWeight)
  
  for (k in 1:length(markPWMs))
  {
    markPWMs[[k]]@name <- paste0(gsub("(\\.bam)?_peaks\\.narrowPeak","",markers[i]),"_","AthLeaf",sep="_",k,sep="_",markPWMs[[k]]@consensus)
  }
  return(markPWMs)
}


getSeqFeature1 <- function(AthFaFile="/home/yanglab/lic/GenomeRef/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
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



getMarkerReadsSignal1 <- function(peakKmerResReads.df,bamFile,peakFile)
{
  library(bamsignals)
  peak <-import(peakFile)
  peak <- peak[seqnames(peak) != "Mt" & seqnames(peak) != "Pt"]
  peakKmerResReads.df$signal <- bamCount(bamFile,peak,verbose=FALSE)
  peakKmerResReads.df$signal <- log2(peakKmerResReads.df$signal +1)
  return(peakKmerResReads.df)
}


cv <- function(peakKmerRes.df, k=10)
{
  require(ranger)
  require(caret)
  folds <- createFolds(peakKmerRes.df$signal, k = k)
  lapply(folds, function(x)
  {
    cv_train <- peakKmerRes.df[x, ]
    cv_test <- peakKmerRes.df[-x, ]
    cv_model <- ranger(signal~., cv_train,num.threads = 6)
    cv_pred <- predict(cv_model, cv_test)
    corRes <- cor.test(cv_pred$prediction,cv_test$signal)
    return(c(corRes$p.value, corRes$estimate, cv_model$r.squared))
  }
  )
}

getImpKmerPerformance <- function(peakFile,bamFile,EnrichImpKmer.df)
{
  mark_pos_KmerRes.df <- getSeqFeature1(peakFile=peakFile, k=6)
  mark_pos_KmerRes.df <- getMarkerReadsSignal1(mark_pos_KmerRes.df, 
                                               bamFile=bamFile,peakFile=peakFile)
  
  mark_pos_KmerRes_EnrichImp.df <- mark_pos_KmerRes.df[,match(rownames(EnrichImpKmer.df),colnames(mark_pos_KmerRes.df)[-2081])]
  
  mark_pos_KmerRes_EnrichImp.df$signal <- mark_pos_KmerRes.df$signal
  
  cc1 <- cv(mark_pos_KmerRes_EnrichImp.df,k = 10)
  return(cc1)
}




AthLeafAllPosNegDf <- vector("list", length = length(markers))
names(AthLeafAllPosNegDf) <- gsub("(\\.bam)?_peaks\\.narrowPeak","",markers)

AthLeafAllLASSOImp <- vector("list", length = length(markers))
names(AthLeafAllLASSOImp) <- gsub("(\\.bam)?_peaks\\.narrowPeak","",markers)

AthLeafAllEnrichmentImp <- vector("list", length = length(markers))
names(AthLeafAllEnrichmentImp) <- gsub("(\\.bam)?_peaks\\.narrowPeak","",markers)

AthLeafAllMismatch_k_mers_list <- vector("list", length = length(markers))
names(AthLeafAllMismatch_k_mers_list) <- gsub("(\\.bam)?_peaks\\.narrowPeak","",markers)

AthLeafmark_pos_KmerRes_EnrichImpCVPer <- vector("list", length = length(markers))
names(AthLeafmark_pos_KmerRes_EnrichImpCVPer) <- gsub("(\\.bam)?_peaks\\.narrowPeak","",markers)


AthLeafAllPWMs<- vector("list", length = length(markers))
names(AthLeafAllPWMs) <- gsub("(\\.bam)?_peaks\\.narrowPeak","",markers)

for (i in 1:length(markers)){
  mark_pos_negDateSet <- getDateSet(peak.postive = markers[i], peak.negtive = markersNeg[i])
  mark_pos_neg_KmerRes.df <- getSeqFeature(peakFile=mark_pos_negDateSet, k=6)
  mark_pos_neg_KmerRes.df <- getMarkerReadsSignal(mark_pos_neg_KmerRes.df, 
                                                  bamFile=bamFiles[i],peakFile=mark_pos_negDateSet)
  AthLeafAllPosNegDf[[i]] <- mark_pos_neg_KmerRes.df
  
  labc <-getActiveVar(peakSignal.df =mark_pos_neg_KmerRes.df)
  AthLeafAllLASSOImp[[i]] <- labc
  aa1 <- getEnrichmentImp(LASSO.impAll = labc,peakFile = markers[i])
  AthLeafAllEnrichmentImp[[i]] <- aa1
  bb1 <- getMismatch_k_mers_list(KmerWeight = aa1$newKmerWeight,enrichment.imp = aa1$EF1_df)
  AthLeafAllMismatch_k_mers_list[[i]] <- bb1
  markPWMs <- getEachMarkPWMs(kmerList = bb1,KmerWeight = aa1$newKmerWeight)
  AthLeafAllPWMs[[i]] <- markPWMs
  cc1 <- getImpKmerPerformance(peakFile=markers[i],bamFile=bamFiles[i],EnrichImpKmer.df=aa1$EF1_df)
  AthLeafmark_pos_KmerRes_EnrichImpCVPer[[i]] <- cc1
  saveRDS(AthLeafAllPosNegDf[[i]],
          file=paste0(gsub("(\\.bam)?_peaks\\.narrowPeak","",markers[i]), "_PosNegDf.RData")) 
  saveRDS(AthLeafAllLASSOImp[[i]],
          file=paste0(gsub("(\\.bam)?_peaks\\.narrowPeak","",markers[i]), "_AthLeafAllLASSOImp.RData")) 
  saveRDS(AthLeafAllEnrichmentImp[[i]],
          file=paste0(gsub("(\\.bam)?_peaks\\.narrowPeak","",markers[i]), "_AthLeafAllEnrichmentImp.RData")) 
  saveRDS(AthLeafAllMismatch_k_mers_list[[i]],
          file=paste0(gsub("(\\.bam)?_peaks\\.narrowPeak","",markers[i]), "_AthLeafAllMismatch_k_mers_list.RData")) 
  
  saveRDS(AthLeafmark_pos_KmerRes_EnrichImpCVPer[[i]],
          file=paste0(gsub("(\\.bam)?_peaks\\.narrowPeak","",markers[i]), "_AthLeafmark_pos_KmerRes_EnrichImpCVPer.RData")) 
  
  write_meme(AthLeafAllPWMs[[i]], 
             file = paste0(gsub("(\\.bam)?_peaks\\.narrowPeak","",markers[i]),"_pos_neg_pwm1.0.meme"))
  print(sprintf("i=%d",i))
}

saveRDS(AthLeafAllPosNegDf,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMs/newET1_newWeight/AthLeafAllPosNegDf.RDS")

saveRDS(AthLeafAllLASSOImp,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMs/newET1_newWeight/AthLeafAllLASSOImp.RDS")
saveRDS(AthLeafAllEnrichmentImp,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMs/newET1_newWeight/AthLeafAllEnrichmentImp.RDS")
saveRDS(AthLeafAllMismatch_k_mers_list,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMs/newET1_newWeight/AthLeafAllMismatch_k_mers_list.RDS")
saveRDS(AthLeafmark_pos_KmerRes_EnrichImpCVPer,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMs/newET1_newWeight/AthLeafmark_pos_KmerRes_EnrichImpCVPer.RDS")

## save the PWMs
saveRDS(AthLeafAllPWMs,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMs/newET1_newWeight/AthLeafAllPWMs.RDS")









