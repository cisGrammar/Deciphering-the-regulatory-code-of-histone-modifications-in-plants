########################################################################################
# Date: 2021 7 24 
# use the fimo files to get the motifs location performance with different binlength & num
# method : Nature Method ; 
# leaf in Arabidopsis thaliana
# pwd : /home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafBinNumDiff
########################################################################################
# binNum = 20
rm(list=ls())
gc()

binNum <- 20

# raw data of each bin
####################################################
setwd("/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafSummitSeq/AthLeafPeakSummit_fimo_thresh1e3")
gffFiles <- list.files(pattern = "gff")
# gffFiles <-gffFiles[seq(1,length(gffFiles),2)]
peakSummit1000 <- list.files(path = "/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafPeakSummit1000_sub",
                             pattern = "bed",full.names = T)
# step_1 :

getMotifLocation <-function(gffFile,peakFile,binNum){
  library(rtracklayer)
  library(GenomicRanges)
  motif.gr <- import.gff3(gffFile)
  motif.gr$motifID <- sapply(strsplit(motif.gr$ID, split = "-"), `[[`, 1)
  peaks.df <- read.table(peakFile,
                         header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  peaks.df <- subset(peaks.df,peaks.df$V1 != "Mt" & peaks.df$V1 != "Pt")
  peaks.df$id <- sprintf("%s:%d-%d", peaks.df$V1, peaks.df$V2, peaks.df$V3)
  peaks.gr <- with(peaks.df, GRanges(seqnames = id, IRanges(start = 1, end = V3 - V2), name = V4))
  bins.gr <- tile(peaks.gr, n = binNum)
  bins.gr <- unlist(bins.gr)
  bins.gr$binID <- 1:binNum
  matPos <- 
    sapply(unique(motif.gr$motifID), function(id){
      motif <- motif.gr[motif.gr$motifID == id]
      sapply(1:binNum, function(i){
        bin <- bins.gr[bins.gr$binID == i]
        sum(countOverlaps(motif, bin, ignore.strand=TRUE))
      })
    })
  storage.mode(matPos) <- "integer"
  return(matPos)
}

############ binNum= 20

AthLeafAllMotifsLocation <- vector("list",length = length(peakSummit1000))
names(AthLeafAllMotifsLocation) <- gsub(".gff","",gffFiles)

for (i in 1:length(gffFiles)){
  AthLeafAllMotifsLocation[[i]] <- getMotifLocation(gffFile = gffFiles[i],
                                                    peakFile = peakSummit1000[i],
                                                    binNum = 20)
  print(i)
  #  saveRDS(AthLeafAllMotifsLocation[[i]],
  #          file=paste0(gsub(".gff","",gffFiles[i]), "_motifLocation.RDS")) 
  gc()
}


## save the AthLeafAllMotifsLocation of bin raw data (bin size 100)
saveRDS(AthLeafAllMotifsLocation,
        file = "/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafBinNumDiff/rawDataBinPlot/AthLeafAllMotifs_peakSummit1000_Location_binNum20RawData.RDS")

############################################## plot locatin #####################
# AthLeafAllMotifsLocation<- readRDS("/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafBinNumDiff/rawDataBinPlot/AthLeafAllMotifs_peakSummit1000_Location_binNum20RawData.RDS")


AthLeafAllMotifsLocation<- readRDS("/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafBinNumDiff/rawDataBinPlot/AthLeafAllMotifs_peakSummit1000_Location_binNum20RawData.RDS")

AthLeafAllMotifsLocation_1 <- AthLeafAllMotifsLocation
for(i in 1:length(AthLeafAllMotifsLocation)){
  AthLeafAllMotifsLocation[[i]] <- t(AthLeafAllMotifsLocation[[i]])
}


binNum <- 20
MarkNum <- sapply(AthLeafAllMotifsLocation,function(x) nrow(x))

MarkName <- names(MarkNum)
for (i in 1:length(AthLeafAllMotifsLocation)){
  AthLeafAllMotifsLocation[[i]] <- data.frame(AthLeafAllMotifsLocation[[i]])
  AthLeafAllMotifsLocation[[i]]$Mark <- MarkName[i]
}



aaa7 <- do.call(rbind,AthLeafAllMotifsLocation)
l1528 <- rownames(aaa7)
l1533 <- aaa7$Mark
names(l1533) <- l1528

split <- factor(l1533, levels=MarkName)

aaa71 <- aaa7[,1:binNum]

aaa8 <-matrix(nrow=nrow(aaa71), ncol = ncol(aaa71))
rownames(aaa8) <-rownames(aaa71)
colnames(aaa8) <-colnames(aaa71)
for (i in 1:nrow(aaa71)) {
  for (j in 1:ncol(aaa71)){
    aaa8[i,j] <- round(aaa71[i,j]/sum(aaa71[i,]),digits = 3)
  }
}

aaa9 <- as.data.frame(aaa8)
aaa9$PWMs <- unlist(lapply(strsplit(rownames(aaa9),"\\."),function(x) x[2]))
markPWMsNum <- as.numeric(table(aaa7$Mark))
total <-  readRDS("/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafAllPWMGC_AT_importance.RDS")
total1 <- merge(aaa9,total,by="PWMs")




library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
cols <- distinctColorPalette(length(AthLeafAllMotifsLocation))
names(cols) <- MarkName
ha = rowAnnotation(df = data.frame(Mark = MarkName),
                   col = list(Mark = cols),
                   width = unit(1, "cm")
)
ha1 = rowAnnotation(df = data.frame(Mark = l1533),
                    col = list(Mark = rep(cols,time=sapply(AthLeafAllMotifsLocation,function(x) nrow(x)))),
                    width = unit(1, "cm"))
mat = as.matrix(aaa8)
colnames(mat) <-c("-1kb","","","","","","","","","Summit","","","","","","","","","","1kb")
rownames(mat) <- rownames(aaa71)
name1 <- c("Depleted","Neutral","Enriched")






Heatmap(mat, name = '     Motif \n enrichment', width = unit(200, "mm"),
        split=split, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        #        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(col=cols,fill=cols))),
        show_row_names = FALSE,
        row_title = NULL,
        column_title = "Bin size 100",
        #        heatmap_legend_param = list(at = c(0.10, 0.2, 0.26), labels = name1)
        right_annotation=ha1) + 
  Heatmap(total1$MotifCatalog, col = c("darkgreen","orange","white"),name = "Motif category", width = unit(5, "mm"))

# bin size 100 binNum 20
# save : /home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafBinNumDiff/rawDataBinPlot/AthLeafAllMotifs_peakSummit1000_Location_binNum20RawData.pdf
# size 10 *7

aaa10 <- as.data.frame(aaa8)
aaa10$Mark <- aaa7$Mark
# plot mean line
dat <- aaa10
colnames(dat)[1:20] <- c("-1kb","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","1kb")
library(reshape2)
dat = melt(dat,variable.name="Bin",value.name = "Ration")
head(dat)
ggplot(data = dat, aes(x = Bin,y = Ration)) + 
  geom_boxplot(fill="grey",width=0.4,notch=F)+
  scale_fill_viridis_c() + facet_wrap(.~Mark, nrow = 3, scales = "free") +
  theme_bw() + theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(0,0.12))
# size : 15 * 7/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPeakSummitLocation/AthLeafBinNumDiff/rawDataBinPlot/AthLeafAllMotifs_peakSummit1000_Location_binNum20RawData_boxPLot.pdf
# save : 
