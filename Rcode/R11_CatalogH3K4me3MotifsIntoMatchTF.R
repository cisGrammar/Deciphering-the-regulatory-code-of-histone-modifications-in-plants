# within AthH3K4me3 peak file
# catalog Motif into match or unMatch TF
# 2021 7 31
###################3
#####################################
# this Rscript compare with /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/AthH3k4me3MotifBedWithinPeak.R

# R
library(rtracklayer)
AthPhastOverlapwithH3K4me3_Genome_motifsBed <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/data/AthPhastOverlapwithH3K4me3_Genome_motifsBed.bedGraph")


############################################3         within AthH3K4me3 peak file

library(universalmotif)
Athh3k4me3PWMs <- read_meme("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/motifs/AthLeafH3K4me3_mis1.meme")
AthH3K4me3MatchTF<- read.table("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/NEWet1.0Mis1/tomtom_H3K4me3_output_files__pwm1.0_overlap5/tomtom.tsv",
                               stringsAsFactors =FALSE,header = TRUE )[,c(1,2,5)]

AthH3K4me3MatchTFNames <- unique(AthH3K4me3MatchTF$Query_ID)
Athh3k4me3UnMatchTFNames <- setdiff(sapply(Athh3k4me3PWMs,function(x) x@name),AthH3K4me3MatchTFNames)


# match TF
# AthRice : 2      total : 12
# AthMaize : 4      total : 16
# AthConservation : 6     total : 21
# AthnoConservation : 13   total : 73


# AthRice


AthRice1 <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/withinPeakMotifBed/AthH3K4me3_AthRiceMotifsGRanges_withinPeak.bed")
AthRice1_match <- AthRice1[which(AthRice1$name %in% AthH3K4me3MatchTFNames)]
AthRice1PhastCon1_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthRice1_match)
AthRice1PhastCon_df1_match <-data.frame("score"= AthRice1PhastCon1_match$score,
                                  "Catalog" = "AthRiceMotifs_match") 

AthRice1_unMatch <- AthRice1[-which(AthRice1$name %in% AthH3K4me3MatchTFNames)]
AthRice1PhastCon1_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthRice1_unMatch)
AthRice1PhastCon_df1_unMatch <-data.frame("score"= AthRice1PhastCon1_unMatch$score,
                                        "Catalog" = "AthRiceMotifs_unMatch") 
length(which(unique(AthRice1$name) %in% AthH3K4me3MatchTFNames))
# match: 2


# AthMaize


AthMaize1 <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/withinPeakMotifBed/AthH3K4me3_AthMaizeMotifsGRanges_withinPeak.bed")
AthMaize1_match <- AthMaize1[which(AthMaize1$name %in% AthH3K4me3MatchTFNames)]
AthMaize1PhastCon1_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthMaize1_match)
AthMaize1PhastCon_df1_match <-data.frame("score"= AthMaize1PhastCon1_match$score,
                                        "Catalog" = "AthMaizeMotifs_match") 

AthMaize1_unMatch <- AthMaize1[-which(AthMaize1$name %in% AthH3K4me3MatchTFNames)]
AthMaize1PhastCon1_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthMaize1_unMatch)
AthMaize1PhastCon_df1_unMatch <-data.frame("score"= AthMaize1PhastCon1_unMatch$score,
                                          "Catalog" = "AthMaizeMotifs_unMatch") 
length(which(unique(AthMaize1$name) %in% AthH3K4me3MatchTFNames))
# match: 4


# Ath conservation


AthConser1 <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/withinPeakMotifBed/AthH3K4me3_conservedMotifsGRanges_withinPeak.bed")
AthConser1_match <- AthConser1[which(AthConser1$name %in% AthH3K4me3MatchTFNames)]
AthConser1PhastCon1_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthConser1_match)
AthConser1PhastCon_df1_match <-data.frame("score"= AthConser1PhastCon1_match$score,
                                        "Catalog" = "AthConserMotifs_match") 

AthConser1_unMatch <- AthConser1[-which(AthConser1$name %in% AthH3K4me3MatchTFNames)]
AthConser1PhastCon1_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthConser1_unMatch)
AthConser1PhastCon_df1_unMatch <-data.frame("score"= AthConser1PhastCon1_unMatch$score,
                                          "Catalog" = "AthConserMotifs_unMatch") 
length(which(unique(AthConser1$name) %in% AthH3K4me3MatchTFNames))
# match: 6



# Ath no-conservation


AthnoConser1 <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/withinPeakMotifBed/AthH3K4me3_noConservedMotifs_withinPeak.bed")
AthnoConser1_match <- AthnoConser1[which(AthnoConser1$name %in% AthH3K4me3MatchTFNames)]
AthnoConser1PhastCon1_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthnoConser1_match)
AthnoConser1PhastCon_df1_match <-data.frame("score"= AthnoConser1PhastCon1_match$score,
                                        "Catalog" = "AthnoConserMotifs_match") 

AthnoConser1_unMatch <- AthnoConser1[-which(AthnoConser1$name %in% AthH3K4me3MatchTFNames)]
AthnoConser1PhastCon1_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthnoConser1_unMatch)
AthnoConser1PhastCon_df1_unMatch <-data.frame("score"= AthnoConser1PhastCon1_unMatch$score,
                                          "Catalog" = "AthnoConserMotifs_unMatch") 
length(which(unique(AthnoConser1$name) %in% AthH3K4me3MatchTFNames))
# match: 13


# Ath random

Athrandom1 <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/withinPeakMotifBed/AthH3K4me3_random_withinPeak.bed")
Ath_PhastCons.bedGraph <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/data/Ath_PhastCons.bedGraph")

Athrandom1PhastCon1 <- subsetByOverlaps(Ath_PhastCons.bedGraph,Athrandom1)

Athrandom1PhastCon_df1 <-data.frame("score"= Athrandom1PhastCon1$score,
                                    "Catalog" = "Background") 


AthAllPWMsCatlogPhastConScore1_TF_withinPeak <- rbind( AthConser1PhastCon_df1_match,AthConser1PhastCon_df1_unMatch,
                                                       AthRice1PhastCon_df1_match,AthRice1PhastCon_df1_unMatch,
                                                       AthMaize1PhastCon_df1_match,AthMaize1PhastCon_df1_unMatch,
                                                       AthnoConser1PhastCon_df1_match,AthnoConser1PhastCon_df1_unMatch,
                                                       Athrandom1PhastCon_df1)
# saveRDS(AthAllPWMsCatlogPhastConScore1_TF_withinPeak,
#          "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/withinPeakMotifBed/boxplot6bp/AthAllPWMsCatlog_TF_WithinPeak_6bp_phastCon.RDS")
        















#################################################3# all AthH3K4me3 catalog all binding site     ##############################333
# catalog with TF

# AthRice

AthRice_genome <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/data/AthH3K4me3_AthRiceMotifsGRanges.bed")
AthRice_genome_match <- AthRice_genome[which(AthRice_genome$name %in% AthH3K4me3MatchTFNames)]
AthRicePhastCon_genome_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthRice_genome_match)
AthRicePhastCon_df_genome_match <-data.frame("score"= AthRicePhastCon_genome_match$score,
                                        "Catalog" = "AthRiceMotifs_match") 

AthRice_genome_unMatch <- AthRice_genome[-which(AthRice_genome$name %in% AthH3K4me3MatchTFNames)]
AthRicePhastCon_genome_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthRice_genome_unMatch)
AthRicePhastCon_df_genome_unMatch <-data.frame("score"= AthRicePhastCon_genome_unMatch$score,
                                               "Catalog" = "AthRiceMotifs_unMatch") 

# AthMaize
AthMaize_genome <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/data/AthH3K4me3_AthMaizeMotifsGRanges.bed")
AthMaize_genome_match <- AthMaize_genome[which(AthMaize_genome$name %in% AthH3K4me3MatchTFNames)]
AthMaizePhastCon_genome_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthMaize_genome_match)
AthMaizePhastCon_df_genome_match <-data.frame("score"= AthMaizePhastCon_genome_match$score,
                                             "Catalog" = "AthMaizeMotifs_match") 

AthMaize_genome_unMatch <- AthMaize_genome[-which(AthMaize_genome$name %in% AthH3K4me3MatchTFNames)]
AthMaizePhastCon_genome_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthMaize_genome_unMatch)
AthMaizePhastCon_df_genome_unMatch <-data.frame("score"= AthMaizePhastCon_genome_unMatch$score,
                                               "Catalog" = "AthMaizeMotifs_unMatch") 



# Ath conservation

AthConser_genome <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/data/AthH3K4me3_conservedMotifsGRanges.bed")
AthConser_genome_match <- AthConser_genome[which(AthConser_genome$name %in% AthH3K4me3MatchTFNames)]
AthConserPhastCon_genome_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthConser_genome_match)
AthConserPhastCon_df_genome_match <-data.frame("score"= AthConserPhastCon_genome_match$score,
                                              "Catalog" = "AthConserMotifs_match") 

AthConser_genome_unMatch <- AthConser_genome[-which(AthConser_genome$name %in% AthH3K4me3MatchTFNames)]
AthConserPhastCon_genome_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthConser_genome_unMatch)
AthConserPhastCon_df_genome_unMatch <-data.frame("score"= AthConserPhastCon_genome_unMatch$score,
                                                "Catalog" = "AthConserMotifs_unMatch")  




# Ath no-conservation

AthnoConser_genome <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/data/AthH3K4me3_noConservedMotifs.bed")
AthnoConser_genome_match <- AthnoConser_genome[which(AthnoConser_genome$name %in% AthH3K4me3MatchTFNames)]
AthnoConserPhastCon_genome_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthnoConser_genome_match)
AthnoConserPhastCon_df_genome_match <-data.frame("score"= AthnoConserPhastCon_genome_match$score,
                                             "Catalog" = "AthnoConserMotifs_match") 

AthnoConser_genome_unMatch <- AthnoConser_genome[-which(AthnoConser_genome$name %in% AthH3K4me3MatchTFNames)]
AthnoConserPhastCon_genome_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthnoConser_genome_unMatch)
AthnoConserPhastCon_df_genome_unMatch <-data.frame("score"= AthnoConserPhastCon_genome_unMatch$score,
                                               "Catalog" = "AthnoConserMotifs_unMatch") 


# Ath random

Athrandom_genome <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/data/Athrandom_Genome.bed")
AthrandomPhastCon_genome <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,Athrandom_genome)

Athrandom1PhastCon_genome_df <-data.frame("score"= AthrandomPhastCon_genome$score,
                                    "Catalog" = "Background") 


AthAllPWMsCatlogPhastConScore1_TF_genome <- rbind(AthConserPhastCon_df_genome_match,AthConserPhastCon_df_genome_unMatch,
                                                  AthRicePhastCon_df_genome_match,AthRicePhastCon_df_genome_unMatch,
                                                 AthMaizePhastCon_df_genome_match,AthMaizePhastCon_df_genome_unMatch,
                                                 AthnoConserPhastCon_df_genome_match,AthnoConserPhastCon_df_genome_unMatch,
                                                Athrandom1PhastCon_genome_df)
# saveRDS(AthAllPWMsCatlogPhastConScore1_TF_genome,
#          "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/Genome/GenomeMotifBed/boxplot6bp/AthAllPWMsCatlog_TF_genome_6bp_phastCon.RDS")















################################################33  intergeneic region
# 

# catalog with TF

# AthRice

AthRice_intergeneic <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/intergenic/intergeneicsPeakMotifBed/AthH3K4me3_AthRiceMotifsGRanges_intergeneic.bed")
AthRice_intergeneic_match <- AthRice_intergeneic[which(AthRice_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthRicePhastCon_intergeneic_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthRice_intergeneic_match)
AthRicePhastCon_df_intergeneic_match <-data.frame("score"= AthRicePhastCon_intergeneic_match$score,
                                             "Catalog" = "AthRiceMotifs_match") 

AthRice_intergeneic_unMatch <- AthRice_intergeneic[-which(AthRice_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthRicePhastCon_intergeneic_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthRice_intergeneic_unMatch)
AthRicePhastCon_df_intergeneic_unMatch <-data.frame("score"= AthRicePhastCon_intergeneic_unMatch$score,
                                               "Catalog" = "AthRiceMotifs_unMatch") 

# AthMaize
AthMaize_intergeneic <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/intergenic/intergeneicsPeakMotifBed/AthH3K4me3_AthMaizeMotifsGRanges_intergeneic.bed")
AthMaize_intergeneic_match <- AthMaize_intergeneic[which(AthMaize_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthMaizePhastCon_intergeneic_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthMaize_intergeneic_match)
AthMaizePhastCon_df_intergeneic_match <-data.frame("score"= AthMaizePhastCon_intergeneic_match$score,
                                              "Catalog" = "AthMaizeMotifs_match") 

AthMaize_intergeneic_unMatch <- AthMaize_intergeneic[-which(AthMaize_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthMaizePhastCon_intergeneic_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthMaize_intergeneic_unMatch)
AthMaizePhastCon_df_intergeneic_unMatch <-data.frame("score"= AthMaizePhastCon_intergeneic_unMatch$score,
                                                "Catalog" = "AthMaizeMotifs_unMatch") 



# Ath conservation

AthConser_intergeneic <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/intergenic/intergeneicsPeakMotifBed/AthH3K4me3_conservedMotifsGRanges_intergeneic.bed")
AthConser_intergeneic_match <- AthConser_intergeneic[which(AthConser_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthConserPhastCon_intergeneic_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthConser_intergeneic_match)
AthConserPhastCon_df_intergeneic_match <-data.frame("score"= AthConserPhastCon_intergeneic_match$score,
                                               "Catalog" = "AthConserMotifs_match") 

AthConser_intergeneic_unMatch <- AthConser_intergeneic[-which(AthConser_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthConserPhastCon_intergeneic_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthConser_intergeneic_unMatch)
AthConserPhastCon_df_intergeneic_unMatch <-data.frame("score"= AthConserPhastCon_intergeneic_unMatch$score,
                                                 "Catalog" = "AthConserMotifs_unMatch")  




# Ath no-conservation

AthnoConser_intergeneic <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/intergenic/intergeneicsPeakMotifBed/AthH3K4me3_noConservedMotifs_intergeneic.bed")
AthnoConser_intergeneic_match <- AthnoConser_intergeneic[which(AthnoConser_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthnoConserPhastCon_intergeneic_match <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthnoConser_intergeneic_match)
AthnoConserPhastCon_df_intergeneic_match <-data.frame("score"= AthnoConserPhastCon_intergeneic_match$score,
                                                 "Catalog" = "AthnoConserMotifs_match") 

AthnoConser_intergeneic_unMatch <- AthnoConser_intergeneic[-which(AthnoConser_intergeneic$name %in% AthH3K4me3MatchTFNames)]
AthnoConserPhastCon_intergeneic_unMatch <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,AthnoConser_intergeneic_unMatch)
AthnoConserPhastCon_df_intergeneic_unMatch <-data.frame("score"= AthnoConserPhastCon_intergeneic_unMatch$score,
                                                   "Catalog" = "AthnoConserMotifs_unMatch") 


# Ath random

Athrandom_intergeneic <- import("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/intergenic/intergeneicsPeakMotifBed/AthH3K4me3_random_intergeneic.bed")
AthrandomPhastCon_intergeneic <- subsetByOverlaps(AthPhastOverlapwithH3K4me3_Genome_motifsBed,Athrandom_intergeneic)

Athrandom1PhastCon_intergeneic_df <-data.frame("score"= AthrandomPhastCon_intergeneic$score,
                                          "Catalog" = "Background") 


AthAllPWMsCatlogPhastConScore1_TF_intergeneic <- rbind(AthConserPhastCon_df_intergeneic_match,AthConserPhastCon_df_intergeneic_unMatch,
                                                  AthRicePhastCon_df_intergeneic_match,AthRicePhastCon_df_intergeneic_unMatch,
                                                  AthMaizePhastCon_df_intergeneic_match,AthMaizePhastCon_df_intergeneic_unMatch,
                                                  AthnoConserPhastCon_df_intergeneic_match,AthnoConserPhastCon_df_intergeneic_unMatch,
                                                  Athrandom1PhastCon_intergeneic_df)
# saveRDS(AthAllPWMsCatlogPhastConScore1_TF_intergeneic,
#          "/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/intergenic/intergeneicsPeakMotifBed/boxplot6bp/AthAllPWMsCatlog_TF_intergeneic_6bp_phastCon.RDS")


############## withinPeak , genome , intergeneic sum

AthAllPWMsCatlogPhastConScore_TF_withinPeak <- readRDS("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/withinPeak/withinPeakMotifBed/boxplot6bp/AthAllPWMsCatlog_TF_WithinPeak_6bp_phastCon.RDS")
AthAllPWMsCatlogPhastConScore_TF_withinPeak$Region <- "WithinPeak"

AthAllPWMsCatlogPhastConScore_TF_genome <- readRDS("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/Genome/GenomeMotifBed/boxplot6bp/AthAllPWMsCatlog_TF_genome_6bp_phastCon.RDS")
AthAllPWMsCatlogPhastConScore_TF_genome$Region <- "Genome"

AthAllPWMsCatlogPhastConScore_TF_intergeneic <- readRDS("/home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/intergenic/intergeneicsPeakMotifBed/boxplot6bp/AthAllPWMsCatlog_TF_intergeneic_6bp_phastCon.RDS")
AthAllPWMsCatlogPhastConScore_TF_intergeneic$Region <- "Intergeneic"

summaryDatas_TF <- rbind(AthAllPWMsCatlogPhastConScore_TF_withinPeak,
                      AthAllPWMsCatlogPhastConScore_TF_genome,
                      AthAllPWMsCatlogPhastConScore_TF_intergeneic)


mean_TF <- aggregate(summaryDatas_TF$score, by=list(summaryDatas_TF$Catalog, summaryDatas_TF$Region), FUN=mean)
len_TF <- aggregate(summaryDatas_TF$score, by=list(summaryDatas_TF$Catalog, summaryDatas_TF$Region), FUN=length)
sd_TF <- aggregate(summaryDatas_TF$score, by=list(summaryDatas_TF$Catalog, summaryDatas_TF$Region), FUN=sd)

mean_TF <- data.frame(mean_TF,len=len_TF$x,sd=sd_TF$x)
colnames(mean_TF) = c( "Catalog","Region", "MeanPhastCon","count","Sd")
mean_TF$Se <- mean_TF$Sd/sqrt(mean_TF$count) ### 计算标准�?



library(RColorBrewer)
library(ggplot2)
mycol= brewer.pal(n = 12, name = "Set3")
ggplot(mean_TF, aes(x=Region, y=MeanPhastCon, fill=Catalog)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.6) +
  scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=MeanPhastCon-Se, ymax=MeanPhastCon +Se),position=position_dodge(.6), width=.2) +
  theme_bw()+
  theme(panel.grid=element_blank()) 

# save : /home/yanglab/Allwork/rework1/Ath_Rice_Mazie_MotifCompare/AthRiceMaizePWMsCompare_yang/l0730boxPlot/CatalogH3K4me3MotifsIntoMatchTF/catalogH3K4me3MotifsIntoMatchTF.pdf
# size : 9 * 5




