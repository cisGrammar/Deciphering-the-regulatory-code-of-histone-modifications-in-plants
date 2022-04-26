# leaf in Arabidopsis thaliana
setwd("/home/yanglab/data/Markers/AthNarrowpeak/leaf/new")   

markers <- list.files(pattern = "bam_peaks.narrowPeak")
bamFiles <- list.files(path = "/home/yanglab/data/Markers/AthNarrowpeak/leaf/new/bam/",
                       pattern = "bam$",full.names = T)

collapseSeq <- function(k = 6,
                        nucleotides = c("A", "T", "C", "G"))
{
  require(gtools)
  require(Biostrings)
  allperms <-
    apply(gtools::permutations(length(nucleotides), k, nucleotides, repeats.allowed = T),
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
  require(ranger)
  require(caret)
  folds <- createFolds(peakKmerResReads.df$signal, k = 10)
  predicted_ten_fold<-c()
  measured_ten_fold<-c()
  for (k in 1:length(folds))
  {
    cv_train <- peakKmerResReads.df[-folds[[k]],  ]
    cv_test <- peakKmerResReads.df[folds[[k]], ]
    cv_model <- ranger(signal~., cv_train)
    cv_pred <- predict(cv_model, cv_test)
    predicted_ten_fold<-c(predicted_ten_fold,cv_pred$predictions)
    measured_ten_fold<-c(measured_ten_fold,cv_test$signal)
  }
  PredictPerformance <- data.frame(predicted = predicted_ten_fold,
                                   observed = measured_ten_fold)
  return(PredictPerformance)
}



cv <- function(peakKmerRes.df, k=10)
{
  require(ranger)
  require(caret)
  folds <- createFolds(peakKmerRes.df$signal, k = k)
  lapply(folds, function(x)
  {
    cv_train <- peakKmerRes.df[-x, ]
    cv_test <- peakKmerRes.df[x, ]
    cv_model <- ranger(signal~., cv_train)
    cv_pred <- predict(cv_model, cv_test)
    corRes <- cor.test(cv_pred$prediction,cv_test$signal)
    return(c(corRes$p.value, corRes$estimate, cv_model$r.squared))
  }
  )
}



AthLeafPerformance_R_PCC <- vector("list", length = length(markers))
names(AthLeafPerformance_R_PCC) <- markers

AthLeafPerformances_obs_pre <- vector("list", length = length(markers))
names(AthLeafPerformances_obs_pre) <- markers

AthLeafAllKmer<- vector("list", length = length(markers))
names(AthLeafAllKmer) <- markers

## get the performance of each mark

for (i in 1:length(bamFiles))
{
  peakKmerResReads.df <- getSeqFeature(peakFile=markers[i], k=6)
  peakKmerResReads.df <- getMarkerReadsSignal(peakKmerResReads.df=peakKmerResReads.df, 
                                              bamFile=bamFiles[i],peakFile = markers[i])
  AthLeafAllKmer[[i]] <- peakKmerResReads.df
  AthLeafPerformance_R_PCC[[i]] <- cv(peakKmerRes.df = peakKmerResReads.df)
  AthLeafPerformances_obs_pre[[i]] <- cv_M_P(peakKmerResReads.df = peakKmerResReads.df)
}


## save the data
saveRDS(AthLeafPerformance_R_PCC,
        file = "/home/yanglab/Allwork/part1/results/Ath/leaf/Rf/AthLeafPerformance_R_PCC.RDS")

saveRDS(AthLeafPerformances_obs_pre,
        file = "/home/yanglab/Allwork/part1/results/Ath/leaf/Rf/AthLeafPerformances_obs_pre.RDS")

saveRDS(AthLeafAllKmer,
        file = "/home/yanglab/Allwork/part1/results/Ath/leaf/Rf/AthLeafAllKmer.RDS")


#######################################################################################
##################################################################################
#############################
# plot the performance of Rsquared and PCC
getmarkersPerformance.Rsquared.melt1 <- function(performance_reads)
{
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- plyr::rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  markersPerformance <-performance_reads
  markersPerformance.tmp <-as.data.frame(t(as.data.frame(markersPerformance)))
  markersPerformance.tmp$marker <- gsub("(\\.bam)?_peaks\\.narrowPeak\\.Fold.*","",rownames(markersPerformance.tmp))
  markersPerformance.tmp$nFold <- gsub(".*(\\.bam)?_peaks\\.narrowPeak\\.","",rownames(markersPerformance.tmp))
  colnames(markersPerformance.tmp) <- c("pvalue", "cor", "Rsquared", "marker","nfold")
  markersPerformance.Rsquared <- data_summary(markersPerformance.tmp[,c(-1,-2)], varname = "Rsquared",groupnames=c("marker"))
  names(markersPerformance.Rsquared) <- c("Marker", "Rsquared", "Rsquared.sd")
  
  markersPerformance.cor <- data_summary(markersPerformance.tmp[,c(-1,-3)], varname = "cor",groupnames=c("marker") )
  markersPerformance.Rsquared$cor <- markersPerformance.cor$cor
  markersPerformance.Rsquared$cor.sd <- markersPerformance.cor$sd
  
  library(reshape2)
  markersPerformance.Rsquared.melt <- melt(markersPerformance.Rsquared, id=c("Marker","Rsquared.sd","cor.sd"))
  markersPerformance.Rsquared1 <- markersPerformance.Rsquared
  names(markersPerformance.Rsquared1) <- c("Marker",  "Rsquared", "sd", "cor", "sd")
  
  
  markersPerformance.Rsquared.melt<-melt(markersPerformance.Rsquared)
  markersPerformance.Rsquared.melt$sd <-0
  markersPerformance.Rsquared.melt[markersPerformance.Rsquared.melt$variable=="Rsquared",]$sd <- markersPerformance.Rsquared.melt[markersPerformance.Rsquared.melt$variable=="Rsquared.sd",]$value
  markersPerformance.Rsquared.melt[markersPerformance.Rsquared.melt$variable=="cor",]$sd <- markersPerformance.Rsquared.melt[markersPerformance.Rsquared.melt$variable=="cor.sd",]$value
  markersPerformance.Rsquared.melt1 <- markersPerformance.Rsquared.melt[markersPerformance.Rsquared.melt$variable=="Rsquared" | markersPerformance.Rsquared.melt$variable=="cor",]
}

######## exclude H2A 2021 7 10
AthLeafPerformance_R_PCC <- AthLeafPerformance_R_PCC[-1]

AthLeafPerformance_R_PCC_plot <- getmarkersPerformance.Rsquared.melt1(AthLeafPerformance_R_PCC)
library(ggplot2)
p<- ggplot(AthLeafPerformance_R_PCC_plot, aes(x=Marker, y= value, fill=variable)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9))

AthLeafPerformance_R_PCCPlot <- p+labs(x = "Histone Marks",y = "Performance")+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(sec.axis = dup_axis(name="Correlation")) +
  scale_fill_discrete(name="Type") + 
  theme(axis.text.x = element_text(angle=90, hjust=0.1, vjust=0.1))

# save : /home/yanglab/Allwork/part1/results/Ath/leaf/Rf/AthLeafPerformance_R_PCCPlot_sub.pdf

##### plot2


## plot the Scatter plot to show the observed signal and predicted signal


# list to data.frame

AthLeafPerformances_obs_pre1 <- lapply(1:length(AthLeafPerformances_obs_pre), function(i) {
  marker <- gsub("(\\.bam)?_peaks\\.narrowPeak","",names(AthLeafPerformances_obs_pre)[i])
  AthLeafPerformances_obs_pre[[i]]$Marker <- marker
  AthLeafPerformances_obs_pre[[i]]
})

res <- data.frame()
for(i in 1:length(AthLeafPerformances_obs_pre1))
  res <- rbind(res, AthLeafPerformances_obs_pre1[[i]])


# plot 
p <- ggplot(data = res, aes(x = observed,y = predicted)) + 
  xlab("Observed histone modification levels") +  ##添加x标签 与上相同
  ylab("Predicted histone modification levels") +
  geom_point() +
  stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
  geom_abline(intercept=0,slope=1) +
  scale_fill_viridis_c() + facet_wrap(.~Marker, nrow = 3, scales = "free") +
  theme_bw()

library(dplyr)
df.cor <- ddply(res, .(Marker), function(val) sprintf("PCC==%.3f (p < 2.2e-16)", cor(val$predicted, val$observed)))

pp <- p+geom_text(data=df.cor, aes(x=0, y=15.9, label=V1), parse=TRUE,hjust = 0,vjust=-0.8,size=4)+
  scale_y_continuous(limits=c(0,20))+
  scale_x_continuous(limits = c(0,20))

#########################################################################################
################################################################


library(dplyr)
df.cor <- ddply(res, .(Marker), function(val) sprintf("PCC==%.3f (p < 2.2e-16)", cor(val$predicted, val$observed)))

pp1 <- p+
  scale_y_continuous(limits=c(0,20))+
  scale_x_continuous(limits = c(0,20))





#############################


# plot AthLeafMarkers
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

AthLeafMarkers <-performance

AthLeafMarkers.tmp <-as.data.frame(t(as.data.frame(AthLeafMarkers)))
AthLeafMarkers.tmp$marker <- gsub("(\\.bam)?_peaks\\.narrowPeak\\.Fold.*","",rownames(AthLeafMarkers.tmp))
AthLeafMarkers.tmp$nFold <- gsub(".*(\\.bam)?_peaks\\.narrowPeak\\.","",rownames(AthLeafMarkers.tmp))
colnames(AthLeafMarkers.tmp) <- c("pvalue", "cor", "Rsquared", "marker","nfold")
AthLeafMarkers.Rsquared <- data_summary(AthLeafMarkers.tmp[,c(-1,-2)], varname = "Rsquared",groupnames=c("marker"))
names(AthLeafMarkers.Rsquared) <- c("Marker", "Rsquared", "Rsquared.sd")

AthLeafMarkers.cor <- data_summary(AthLeafMarkers.tmp[,c(-1,-3)], varname = "cor",groupnames=c("marker") )
AthLeafMarkers.Rsquared$cor <- AthLeafMarkers.cor$cor
AthLeafMarkers.Rsquared$cor.sd <- AthLeafMarkers.cor$sd

AthLeafMarkers.Rsquared.melt <- melt(AthLeafMarkers.Rsquared, id=c("Marker","Rsquared.sd","cor.sd"))
AthLeafMarkers.Rsquared1 <- AthLeafMarkers.Rsquared
names(AthLeafMarkers.Rsquared1) <- c("Marker",  "Rsquared", "sd", "cor", "sd")


melt(AthLeafMarkers.Rsquared)
AthLeafMarkers.Rsquared.melt<-melt(AthLeafMarkers.Rsquared)
AthLeafMarkers.Rsquared.melt$sd <-0
AthLeafMarkers.Rsquared.melt[AthLeafMarkers.Rsquared.melt$variable=="Rsquared",]$sd <- AthLeafMarkers.Rsquared.melt[AthLeafMarkers.Rsquared.melt$variable=="Rsquared.sd",]$value
AthLeafMarkers.Rsquared.melt[AthLeafMarkers.Rsquared.melt$variable=="cor",]$sd <- AthLeafMarkers.Rsquared.melt[AthLeafMarkers.Rsquared.melt$variable=="cor.sd",]$value
AthLeafMarkers.Rsquared.melt1 <- AthLeafMarkers.Rsquared.melt[AthLeafMarkers.Rsquared.melt$variable=="Rsquared" | AthLeafMarkers.Rsquared.melt$variable=="cor",]


p<- ggplot(AthLeafMarkers.Rsquared.melt1, aes(x=Marker, y= value, fill=variable)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9))

p+labs(title="Sequence Predicts Performance of Histone Markers", x = "Histone Markers",y = "Performance")+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(sec.axis = dup_axis(name="Correlation")) +
  scale_fill_discrete(name="Type") + 
  theme(axis.text.x = element_text(angle=90, hjust=0.1, vjust=0.1))








### to get the observed and predicted signal
cv_M_P <- function(peakKmerResReads.df, k=10)
{
  require(ranger)
  require(caret)
  folds <- createFolds(peakKmerResReads.df$signal, k = 10)
  predicted_ten_fold<-c()
  measured_ten_fold<-c()
  for (k in 1:length(folds))
  {
    cv_train <- peakKmerResReads.df[-folds[[k]],  ]
    cv_test <- peakKmerResReads.df[folds[[k]], ]
    cv_model <- ranger(signal~., cv_train)
    cv_pred <- predict(cv_model, cv_test)
    predicted_ten_fold<-c(predicted_ten_fold,cv_pred$predictions)
    measured_ten_fold<-c(measured_ten_fold,cv_test$signal)
  }
  PredictPerformance <- data.frame(predicted = predicted_ten_fold,
                                   observed = measured_ten_fold)
  return(PredictPerformance)
}






