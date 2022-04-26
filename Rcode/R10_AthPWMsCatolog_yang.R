# reCatalog Ath H3K4me3 motif
# 2021 7 30



getDFunique <- function(df){
  l521 <- df
  colnames(l521) <- c("id1","id2")
  l521$id1_id2 <- paste0(l521$id1,"_",l521$id2)
  l521$id2_id1 <- paste0(l521$id2,"_",l521$id1)
  i <- 1
  index1 <- length(which(l521$id1_id2 %in% l521$id2_id1))
  while(index1 != 0) 
  {
    index2 <- which(l521$id2_id1 %in% l521$id1_id2[i])
    if (length(index2) != 0)
    {
      l521 <- l521[-index2,]
    }
    i <- i+1
    index1 <- length(which(l521$id1_id2 %in% l521$id2_id1))
  }
  colnames(l521) <- colnames(df)
  return(l521[,c(1,2)]) 
}





Ath_Rice <- read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/dubleCompare/tomtom_Ath_Rice_output_files_pwm1.0_overlap5/tomtom.tsv",
                       stringsAsFactors =FALSE,header = TRUE )[,c(1,2,5)]
Rice_Ath <- read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/dubleCompare/tomtom_Rice_Ath_output_files_pwm1.0_overlap5/tomtom.tsv",
                       stringsAsFactors =FALSE,header = TRUE )[,c(1,2,5)]
Ath_Rice_Merge <- data.frame("Query" = c(Ath_Rice$Query_ID,Rice_Ath$Target_ID),
                             "Target" = c(Ath_Rice$Target_ID,Rice_Ath$Query_ID))
Ath_Rice_Merge_unique <- unique(Ath_Rice_Merge)

# Ath & Maize
Ath_Maize <- read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/dubleCompare/tomtom_Ath_Maize_output_files_pwm1.0_overlap5/tomtom.tsv",
                        stringsAsFactors =FALSE,header = TRUE )[,c(1,2,5)]
Maize_Ath  <- read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/dubleCompare/tomtom_Maize_Ath_output_files_pwm1.0_overlap5/tomtom.tsv",
                         stringsAsFactors =FALSE,header = TRUE )[,c(1,2,5)]
Ath_Maize_Merge <- data.frame("Query" = c(Ath_Maize$Query_ID,Maize_Ath$Target_ID),
                              "Target" = c(Ath_Maize$Target_ID,Maize_Ath$Query_ID))
Ath_Maize_Merge_unique <- unique(Ath_Maize_Merge)


# Rice & Maize
Rice_Maize <- read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/dubleCompare/tomtom_Rice_Maize_output_files_pwm1.0_overlap5/tomtom.tsv",
                         stringsAsFactors =FALSE,header = TRUE )[,c(1,2,5)]
Maize_Rice <- read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/dubleCompare/tomtom_Maize_Rice_output_files_pwm1.0_overlap5/tomtom.tsv",
                         stringsAsFactors =FALSE,header = TRUE )[,c(1,2,5)]
Rice_Maize_Merge <- data.frame("Query" = c(Rice_Maize$Query_ID,Maize_Rice$Target_ID),
                               "Target" = c(Rice_Maize$Target_ID,Maize_Rice$Query_ID),stringsAsFactors = FALSE)

Rice_Maize_Merge_unique <- unique(Rice_Maize_Merge)

############################ catalog   


##### 2021 7 30
# conserved
conservasion1 <- intersect(as.character(Ath_Maize_Merge_unique$Query),as.character(Ath_Rice_Merge_unique$Query))
RiceMazie_toAth_Rice <- as.character(Ath_Rice_Merge_unique$Query)[as.numeric(na.omit(match(as.character(Rice_Maize_Merge_unique$Query),
                                                                                           as.character( Ath_Rice_Merge_unique$Target))))]
MaizeMazie_toAth_Mazie <- as.character(Ath_Maize_Merge_unique$Query)[as.numeric(na.omit(match(as.character(Rice_Maize_Merge_unique$Target),
                                                                                              as.character(Ath_Maize_Merge_unique$Target))))]
conservedMotifs <- unique(c(conservasion1,RiceMazie_toAth_Rice,MaizeMazie_toAth_Mazie))

AthRice_PWMs_all <- unique(as.character(Ath_Rice_Merge_unique$Query)) # 31
AthMaize_PWMs_all <- unique(as.character(Ath_Maize_Merge_unique$Query)) # 34

AthRice_Motifs <- AthRice_PWMs_all[-which(AthRice_PWMs_all %in% conservedMotifs)] # 12
AthMaize_Motifs <- AthMaize_PWMs_all[-which(AthMaize_PWMs_all %in% conservedMotifs)] # 16


similiarMotifs <- unique(c(as.character(Ath_Rice_Merge_unique$Query),
                           as.character(Ath_Maize_Merge_unique$Query))) # 49


noConservedMotifs <- unlist(lapply(AthLeafH3K4me3Motifs, function(x) x@name))[-which(unlist(lapply(AthLeafH3K4me3Motifs, function(x) x@name)) %in% similiarMotifs)]
# 73


library(universalmotif)
library(rtracklayer)
AthLeafH3K4me3Motifs <- read_meme("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/AthLeafH3K4me3_mis1.meme") 
AthLeafH3K4me3GenomeBS <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/data/H3K4me3_Genome_motifs.bed")
AthLeafH3K4me3GenomeBS_1 <- AthLeafH3K4me3GenomeBS[seqnames(AthLeafH3K4me3GenomeBS) %in% as.character(seq(1:5))]

seqlevels(AthLeafH3K4me3GenomeBS_1) <- paste0("chr",seqlevels(AthLeafH3K4me3GenomeBS_1))

start(AthLeafH3K4me3GenomeBS_1) <- start(AthLeafH3K4me3GenomeBS_1)-1
strand(AthLeafH3K4me3GenomeBS_1) <- "*"
# AthLeafH3K4me3GenomeBS_1$score <- 0


### conserved motif (sample all)
conservedMotifsGRangesList  <- vector("list",length(conservedMotifs))
for (i in 1:length(conservedMotifs))
{
  conservedMotifsGRangesList[[i]] <- AthLeafH3K4me3GenomeBS_1[which(AthLeafH3K4me3GenomeBS_1$name == conservedMotifs[i])]
}
conservedMotifsGRanges <- sort.GenomicRanges(unlist(GRangesList(conservedMotifsGRangesList)))
conservedMotifsGRanges <-granges(conservedMotifsGRanges)
# seqlevels(conservedMotifsGRanges) <- paste0("chr",seqlevels(conservedMotifsGRanges))
export(conservedMotifsGRanges,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_conservedMotifsGRanges.bed")

### AthRice motif (sample all)
AthRiceMotifsGRangesList  <- vector("list",length(AthRice_Motifs))
for (i in 1:length(AthRice_Motifs))
{
  AthRiceMotifsGRangesList[[i]] <- AthLeafH3K4me3GenomeBS_1[which(AthLeafH3K4me3GenomeBS_1$name == AthRice_Motifs[i])]
}
AthRiceMotifsGRanges <- sort.GenomicRanges(unlist(GRangesList(AthRiceMotifsGRangesList)))
AthRiceMotifsGRanges <-granges(AthRiceMotifsGRanges)
# seqlevels(AthRiceMotifsGRanges) <- paste0("chr",seqlevels(AthRiceMotifsGRanges))
export(AthRiceMotifsGRanges,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthRiceMotifsGRanges.bed")


### RiceMaize motif (sample all)
AthMaizeMotifsGRangesList  <- vector("list",length(AthMaize_Motifs))
for (i in 1:length(AthMaize_Motifs))
{
  AthMaizeMotifsGRangesList[[i]] <- AthLeafH3K4me3GenomeBS_1[which(AthLeafH3K4me3GenomeBS_1$name == AthMaize_Motifs[i])]
}
AthMaizeMotifsGRanges <- sort.GenomicRanges(unlist(GRangesList(AthMaizeMotifsGRangesList)))
AthMaizeMotifsGRanges <-granges(AthMaizeMotifsGRanges)
# seqlevels(AthMaizeMotifsGRanges) <- paste0("chr",seqlevels(AthMaizeMotifsGRanges))
export(AthMaizeMotifsGRanges,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthMaizeMotifsGRanges.bed")



noConservedMotifsList  <- vector("list",length(noConservedMotifs))
for (i in 1:length(noConservedMotifs))
{
  noConservedMotifsList[[i]] <- AthLeafH3K4me3GenomeBS_1[which(AthLeafH3K4me3GenomeBS_1$name == noConservedMotifs[i])]
}
AthnoConservedMotifs <- sort.GenomicRanges(unlist(GRangesList(noConservedMotifsList)))
AthnoConservedMotifs <-granges(AthnoConservedMotifs)
# seqlevels(AthnoConservedMotifs) <- paste0("chr",seqlevels(AthnoConservedMotifs))
export(AthnoConservedMotifs,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_noConservedMotifs.bed")





###########################################
# bwtool
# conserved motif 21

##########################################################3333##### phastCon
# AthMaize
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthMaizeMotifsGRanges.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/data/Ath_PhastCons.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/AthH3K4me3_AthMaizeMotifsGRanges_mean_signal.txt
# AthRice
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthRiceMotifsGRanges.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/data/Ath_PhastCons.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/AthH3K4me3_AthRiceMotifsGRanges_mean_signal.txt
# conservation
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_conservedMotifsGRanges.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/data/Ath_PhastCons.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/AthH3K4me3_conservedMotifsGRanges_mean_signal.txt




############################ plot phastCon

conservedMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/AthH3K4me3_conservedMotifsGRanges_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)

noConservedMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/newClassification_0709/0728_withoutStrand_allSample_replot/AthH3K4me3_noConservedMotifsGRanges_withoutStrand_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)
# Ath & Rice
AthRiceMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/AthH3K4me3_AthRiceMotifsGRanges_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)

# Rice Maize
AthMaizeMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/AthH3K4me3_AthMaizeMotifsGRanges_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)


AthRandom_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/newClassification_0709/0728_withoutStrand_allSample_replot/Athrandom_withoutStrand_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)

AthLeafH3K4me3PhastCons <- data.frame("Pos" = conservedMotifs_signal$V1,
                                      "conservedMotifs" = conservedMotifs_signal$V2,
                                      "Non-ConservedMotifs" = noConservedMotifs_signal$V2 ,
                                      "AthRiceMotifs" = AthRiceMotifs_signal$V2,
                                      "AthMaizeMotifs" = AthMaizeMotifs_signal$V2,
                                      "Background" = AthRandom_signal$V2)

AthLeafH3K4me3PhastCons <-as.data.frame(apply(AthLeafH3K4me3PhastCons,2,as.numeric))
AthLeafH3K4me3PhastConsMelt <- melt(AthLeafH3K4me3PhastCons, id = "Pos")


ggplot(data = AthLeafH3K4me3PhastConsMelt, aes(x = Pos, y = value,color = variable)) +
  geom_line(size = 0.7) + 
  theme_classic() +
  theme(legend.position = 'top') + 
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + 
  labs(x = '', y = 'PhastCon score') +
  geom_vline(aes(xintercept = -3), color="black",
             linetype="dashed", size=0.4)+
  geom_vline(aes(xintercept = 3), color="black",
             linetype="dashed", size=0.4)
# save : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/H3K4me3_RicePhastCon_Range500_conser21.pdf
# size : 9 * 6

AthLeafH3K4me3PhastCons1 <- AthLeafH3K4me3PhastCons[300:700,]

AthLeafH3K4me3PhastCons1 <-as.data.frame(apply(AthLeafH3K4me3PhastCons1,2,as.numeric))
AthLeafH3K4me3PhastConsMelt1 <- melt(AthLeafH3K4me3PhastCons1, id = "Pos")

ggplot(data = AthLeafH3K4me3PhastConsMelt1, aes(x = Pos, y = value,color = variable)) +
  geom_line(size = 0.7) + 
  theme_classic() +
  theme(legend.position = 'top') + 
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + 
  labs(x = '', y = 'PhastCon score') +
  geom_vline(aes(xintercept = -3), color="black",
             linetype="dashed", size=0.4)+
  geom_vline(aes(xintercept = 3), color="black",
             linetype="dashed", size=0.4)

# save : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phastCon/H3K4me3_RicePhastCon_Range200_conser21.pdf
# size : 9 * 6














############## phylop
# AthMaize
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthMaizeMotifsGRanges.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/phylop/Ath_phylop_noMtPt_chr.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phylop/AthH3K4me3_AthMaizeMotifsGRanges_mean_signal.txt
# AthRice
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthRiceMotifsGRanges.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/phylop/Ath_phylop_noMtPt_chr.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phylop/AthH3K4me3_AthRiceMotifsGRanges_mean_signal.txt
# conservation
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_conservedMotifsGRanges.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/phylop/Ath_phylop_noMtPt_chr.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/phylop/AthH3K4me3_conservedMotifsGRanges_mean_signal.txt



























############### interneics region
# pwd : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/newClassification_intergentic0726/intergenicsRange_lzh
#######################################################
# get no intergentic from Ath GENOME 
# 2021 7 30
##################################

#使用bedtools获得随机序列2000条
## ####################### R 
library(rtracklayer)
AthGenomeIntergenic <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/newClassification_intergentic0726/intergenicsRange_lzh/intergeneic.bed")
seqlevels(AthGenomeIntergenic) <- paste0("chr",seqlevels(AthGenomeIntergenic))

# AthRice
AthRice <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthRiceMotifsGRanges.bed")
AthRice1 <- subsetByOverlaps(AthRice,AthGenomeIntergenic,minoverlap = 6)
export(AthRice1,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/AthH3K4me3_AthRiceMotifsGRanges_interneic.bed")


# AthMaize
AthMaize <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_AthMaizeMotifsGRanges.bed")
AthMaize1 <- subsetByOverlaps(AthMaize,AthGenomeIntergenic,minoverlap = 6)
export(AthMaize1,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/AthH3K4me3_AthMaizeMotifsGRanges_interneic.bed")


# Ath conservation
AthConser <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthH3K4me3_conservedMotifsGRanges.bed")
AthConser1 <- subsetByOverlaps(AthConser,AthGenomeIntergenic,minoverlap = 6)
export(AthConser1,"/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/AthH3K4me3_conservedMotifsGRanges_interneic.bed")

# AthMaize
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/AthH3K4me3_AthMaizeMotifsGRanges_interneic.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/data/Ath_PhastCons.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/AthH3K4me3_AthMaizeMotifsGRanges_mean_signal.txt
# AthRice
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/AthH3K4me3_AthRiceMotifsGRanges_interneic.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/data/Ath_PhastCons.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/AthH3K4me3_AthRiceMotifsGRanges_mean_signal.txt
# conservation
bwtool agg 500:500  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/AthH3K4me3_conservedMotifsGRanges_interneic.bed  /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/data/Ath_PhastCons.bedGraph.bw /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/AthH3K4me3_conservedMotifsGRanges_mean_signal.txt

############################ plot phastCon

conservedMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/AthH3K4me3_conservedMotifsGRanges_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)

noConservedMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/newClassification_intergentic0726/intergenicsRange_lzh/AthH3K4me3_noConservedMotifsGRanges_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)
# Ath & Rice
AthRiceMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/AthH3K4me3_AthRiceMotifsGRanges_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)

# Ath Maize
AthMaizeMotifs_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/AthH3K4me3_AthMaizeMotifsGRanges_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)


AthRandom_signal <-read.table("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/newClassification_intergentic0726/intergenicsRange_lzh/Athrandom_internics_mean_signal.txt",header=F,sep="\t",stringsAsFactors= F)

AthLeafH3K4me3PhastCons <- data.frame("Pos" = conservedMotifs_signal$V1,
                                      "conservedMotifs" = conservedMotifs_signal$V2,
                                      "Non-ConservedMotifs" = noConservedMotifs_signal$V2 ,
                                      "AthRiceMotifs" = AthRiceMotifs_signal$V2,
                                      "AthMaizeMotifs" = AthMaizeMotifs_signal$V2,
                                      "Background" = AthRandom_signal$V2)

AthLeafH3K4me3PhastCons <-as.data.frame(apply(AthLeafH3K4me3PhastCons,2,as.numeric))
AthLeafH3K4me3PhastConsMelt <- melt(AthLeafH3K4me3PhastCons, id = "Pos")


ggplot(data = AthLeafH3K4me3PhastConsMelt, aes(x = Pos, y = value,color = variable)) +
  geom_line(size = 0.7) + 
  theme_classic() +
  theme(legend.position = 'top') + 
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + 
  labs(x = '', y = 'PhastCon score') +
  geom_vline(aes(xintercept = -3), color="black",
             linetype="dashed", size=0.4)+
  geom_vline(aes(xintercept = 3), color="black",
             linetype="dashed", size=0.4)
# save : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/H3K4me3_RicePhastCon_Range500_conser21_intergenic.pdf
# size : 9 * 6

AthLeafH3K4me3PhastCons1 <- AthLeafH3K4me3PhastCons[300:700,]

AthLeafH3K4me3PhastCons1 <-as.data.frame(apply(AthLeafH3K4me3PhastCons1,2,as.numeric))
AthLeafH3K4me3PhastConsMelt1 <- melt(AthLeafH3K4me3PhastCons1, id = "Pos")

ggplot(data = AthLeafH3K4me3PhastConsMelt1, aes(x = Pos, y = value,color = variable)) +
  geom_line(size = 0.7) + 
  theme_classic() +
  theme(legend.position = 'top') + 
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + 
  labs(x = '', y = 'PhastCon score') +
  geom_vline(aes(xintercept = -3), color="black",
             linetype="dashed", size=0.4)+
  geom_vline(aes(xintercept = 3), color="black",
             linetype="dashed", size=0.4)

# save : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/H3K4me3_RicePhastCon_Range200_conser21_intergenic.pdf
# size : 9 * 6

AthLeafH3K4me3PhastCons2 <- AthLeafH3K4me3PhastCons[400:600,]

AthLeafH3K4me3PhastCons2 <-as.data.frame(apply(AthLeafH3K4me3PhastCons2,2,as.numeric))
AthLeafH3K4me3PhastConsMelt2 <- melt(AthLeafH3K4me3PhastCons2, id = "Pos")

ggplot(data = AthLeafH3K4me3PhastConsMelt2, aes(x = Pos, y = value,color = variable)) +
  geom_line(size = 0.7) + 
  theme_classic() +
  theme(legend.position = 'top') + 
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + 
  labs(x = '', y = 'PhastCon score') +
  geom_vline(aes(xintercept = -3), color="black",
             linetype="dashed", size=0.4)+
  geom_vline(aes(xintercept = 3), color="black",
             linetype="dashed", size=0.4)

# save : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/intergenic_lzh/phastCon/H3K4me3_RicePhastCon_Range100_conser21_intergenic.pdf
# size : 9 * 6









































### plot number of different motifs
# 2021 7 21 
length(unique(Rice_MaizeH3K4me3TOM$Query_ID))
# 24
# Ath : conservation   : 21
#       noConservation : 73
#       AthRice        : 12
#       AthMaize       : 16
################################ pie plot 
ad = data.frame(type=c("Conserved","Non-coserved","Arabidopsis_Rice","Arabidopsis_Mazie"),
                value=c(21,73,12,16))
ad$fraction = ad$value / sum(ad$value)
ad$ymax = cumsum(ad$fraction)
ad$ymin = c(0, head(ad$ymax, n = -1))

label_value <- paste('(', round(ad$value/sum(ad$value) * 100, 1), '%)', sep = '')
Mylabel <- paste(ad$type, label_value, sep = '')

Mylabel1 <- c(Mylabel[4],Mylabel[3],Mylabel[1],Mylabel[2])
library(ggplot2)
ggplot(data = ad, aes(fill = type, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(colour = "grey30") +
  coord_polar(theta = "y") +
  labs(x = "", y = "", title = "") + 
  xlim(c(1.5, 4)) +
  theme_bw() +
  theme(panel.grid=element_blank()) + ## 去掉白色外框
  theme(axis.text=element_blank()) + ## 把图旁边的标签去掉
  theme(axis.ticks=element_blank()) + ## 去掉左上角的坐标刻度线
  theme(panel.border=element_blank()) + ## 去掉最外层的正方形边框
  geom_text(aes(x = 3.5, y = ((ymin+ymax)/2), label = value)) +
  scale_fill_discrete(labels = Mylabel1)

# save : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/catalogBed/AthLeafPWMconservationCatalog_pie_yang.pdf
# size :  7.5 * 6.5


# Ath : 122   Rice : 82   Maize : 91
# AthRice : 33   AthMazie : 37  RiceMaize : 24
# AthRiceMaize : 21
library(venneuler)
vd <- euler(c(Arabidopsis = 122, Rice = 82, Mazie = 91,
              "Arabidopsis&Rice" = 33, "Arabidopsis&Mazie" = 37, "Rice&Mazie" = 24,
              "Arabidopsis&Rice&Mazie" = 16))
plot(vd,
     fills = list(fill =c(colors()[616], colors()[38], colors()[468]), alpha = 0.6),
     labels = list(col = "red", font = 4), 
     edges = FALSE,
     quantities = TRUE)


library(VennDiagram)
venn.plot <- draw.triple.venn(
  area1 = 122,
  area2 = 82,
  area3 = 91,
  n12 = 33,
  n23 = 24,
  n13 = 37,
  n123 = 21,
  category = c("Arabidopsis", "Rice", "Mazie"),
  fill =c(colors()[616], colors()[38], colors()[468]),
  lty = "blank")
grid.draw(venn.plot)


venn.plot <- draw.triple.venn(area1 = 122, area2 =82, area3 = 91, 
                              n12 = 33, n23 = 24, n13 = 37, n123 = 21, 
                              category = c("Arabidopsis", "Rice", "Mazie"), 
                              fill = c(colors()[616], colors()[38], colors()[468]),
                              cex = 1,cat.cex = 1.5,
                              ext.line.lty = "dashed" )

# SAVE : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/RicePhylop/data/CatalogBed/Ath_Rice_Maize_H3K4me3PWMs_VEN_yang_6_6.pdf
# size : 6 * 6






