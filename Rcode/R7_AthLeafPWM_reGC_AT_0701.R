########################################
# pwd : /home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafPWM_reGC_AT.R
# 2021 6 29
##################################

library(stringr)
library(universalmotif)

allAthLeafMotifs <- read_meme("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/NEWet1.0Mis1/AthLeafAllPWMs_family.meme")
Ath_TF_binding_motifs <- read_meme("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/NEWet1.0Mis1/Ath_TF_binding_motifs.meme")
### motif family 
for (i in 1:length(allAthLeafMotifs)) allAthLeafMotifs[[i]]@family <- as.character(unlist(strsplit(allAthLeafMotifs[[i]]@name, split = "_")))[1]

for (i in 1:length(Ath_TF_binding_motifs)) 
{
  Ath_TF_binding_motifs[[i]]@multifreq$GC <- (sum(Ath_TF_binding_motifs[[i]]@motif[rownames(Ath_TF_binding_motifs[[i]]@motif) == "G"] + Ath_TF_binding_motifs[[i]]@motif[rownames(Ath_TF_binding_motifs[[i]]@motif) == "C"]))/ncol(Ath_TF_binding_motifs[[i]])
  Ath_TF_binding_motifs[[i]]@multifreq$AT <- (sum(Ath_TF_binding_motifs[[i]]@motif[rownames(Ath_TF_binding_motifs[[i]]@motif) == "A"] + Ath_TF_binding_motifs[[i]]@motif[rownames(Ath_TF_binding_motifs[[i]]@motif) == "T"]))/ncol(Ath_TF_binding_motifs[[i]])
}
#hist(unlist(lapply(Ath_TF_binding_motifs,function(x) x@multifreq$GC)),breaks=20)
Ath_TF_binding_motifs_GC <- unlist(lapply(Ath_TF_binding_motifs,function(x) x@multifreq$GC))
#Ath_TF_binding_motifs_GC <- round((sum(unlist(lapply(Ath_TF_binding_motifs,function(x) x@multifreq$GC)))) / length(Ath_TF_binding_motifs),4)
#mean GC: 0.4555
#hist(unlist(lapply(Ath_TF_binding_motifs,function(x) x@multifreq$AT)),breaks=20)
Ath_TF_binding_motifs_AT <- unlist(lapply(Ath_TF_binding_motifs,function(x) x@multifreq$AT))
#Ath_TF_binding_motifs_AT <- round((sum(unlist(lapply(Ath_TF_binding_motifs,function(x) x@multifreq$AT)))) / length(Ath_TF_binding_motifs),4)
#mean AT 0.5445

################ AthLeaf PWMS

for (i in 1:length(allAthLeafMotifs)) 
{
  allAthLeafMotifs[[i]]@multifreq$GC <- round((sum(allAthLeafMotifs[[i]]@motif[rownames(allAthLeafMotifs[[i]]@motif) == "G"] + allAthLeafMotifs[[i]]@motif[rownames(allAthLeafMotifs[[i]]@motif) == "C"]))/ncol(allAthLeafMotifs[[i]]),4)
  allAthLeafMotifs[[i]]@multifreq$AT <- round((sum(allAthLeafMotifs[[i]]@motif[rownames(allAthLeafMotifs[[i]]@motif) == "A"] + allAthLeafMotifs[[i]]@motif[rownames(allAthLeafMotifs[[i]]@motif) == "T"]))/ncol(allAthLeafMotifs[[i]]),4)
}
allAthLeafMotifs_GC <- unlist(lapply(allAthLeafMotifs,function(x) x@multifreq$GC))
# 0.5325437
allAthLeafMotifs_AT <- unlist(lapply(allAthLeafMotifs,function(x) x@multifreq$AT))
# 0.4674562
GCIndexs <- which(allAthLeafMotifs_GC > quantile(Ath_TF_binding_motifs_GC)[4])
# length(GCIndex) = 646
ATIndexs <- which(allAthLeafMotifs_GC < quantile(Ath_TF_binding_motifs_GC)[2])
# length(ATIndex) = 199
motifCal <- rep("Neutral motifs",length(allAthLeafMotifs))
motifCal[GCIndexs] <- "GC-rich motif"
motifCal[ATIndexs] <- "AT-rich motif"

PWMsNames <- unlist(lapply(allAthLeafMotifs,function(x) x@name))
PWMsMark <- unlist(lapply(allAthLeafMotifs,function(x) x@family))
df1 <- data.frame("PWMs" = PWMsNames,
                  "Mark" = PWMsMark,
                  "MotifCatalog"= motifCal)
allPWMsCatalog <-tapply(df1$Mark, df1$MotifCatalog, table)
PWMsSummary <- data.frame("allMOtif" = as.numeric(table(df1$Mark)),
                          "GC-rich motif" = as.numeric(allPWMsCatalog$`GC-rich motif`),
                          "AT-rich motif" = as.numeric(allPWMsCatalog$`AT-rich motif`),
                          "Neutral motif" = as.numeric(allPWMsCatalog$`Neutral motifs`))

rownames(PWMsSummary) <- unique(as.character(unlist(lapply(allAthLeafMotifs, function(x) x@family))))

#write.table(PWMsSummary,file = "/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafMotifAT_GC_Num.txt",sep = "\t",quote = F,row.names = TRUE,col.names = TRUE)


################################ plot ######################################

AthLeafMotifGC_AT_Num <-barplot(height =t(as.matrix(PWMsSummary[,-1])),  # 绘图数据（矩阵）
                                col = c('orange', 'steelblue','grey'),  # 填充颜色
                                border = '#ffffff',   # 轮廓颜色
                                axisnames = FALSE,
                                ylab = 'Number of motifs',  # Y轴名称
                                horiz = FALSE,  # 是否为水平放置
                                ylim = c(0, 150), # Y轴取值范围
                                legend.text = c('G+C rich motif', 'A+T rich motif', 'Neutral motifs'),  # 图例文本
                                args.legend=c(x=15.5,y=150,cex=0.8),cex.axis = 0.8)
text(AthLeafMotifGC_AT_Num, par("usr")[3], labels = rownames(PWMsSummary), srt = 60, adj = c(1.1,1.1), xpd = TRUE, cex=.8)


############ plot %
data1 <- data.frame("Mark" = rep(rownames(PWMsSummary),times=3),
                    "Catalog" = rep(c('G+C rich motif', 'A+T rich motif', 'Neutral motifs'),each=13),
                    "Value" = c(PWMsSummary$GC.rich.motif,PWMsSummary$AT.rich.motif,PWMsSummary$Neutral.motif))

data1_sub <- subset(data1,Mark != "H2A")
library(ggplot2)
ggplot(data1_sub, aes(fill=Catalog, y=Value, x=Mark)) +
       geom_bar(position="fill", stat="identity") + 
#  scale_fill_viridis(discrete = T,option = "H") +
       theme(axis.text.x = element_text(angle=90, hjust=0.1, vjust=0.1))

# save : /home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafMotifAT_GC_replot_1_sub.pdf

##############################################################
# add model importance
AthLeafAllPWMsImporatntPred<-readRDS("/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafTFfamily/AthLeafAllPWMImportance.RDS")

AthLeafAllPWMsImporatntPred1 <- do.call(rbind,AthLeafAllPWMsImporatntPred)
total <- merge(df1,AthLeafAllPWMsImporatntPred1,by="PWMs")

# saveRDS(total,"/home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafAllPWMGC_AT_importance.RDS")

######### plot weight

ggplot(total,aes(weight,..density..))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))

ggplot(total,aes(weight,..density.., color=MotifCatalog))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))


ggplot(total,aes(weight1,..density..))  +
  geom_boxplot(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) 


############## SOFTWARE 
softmax2 <- function(x, scaleFactor=0.6){
  tmp <- log10(x) / scaleFactor
  return(exp(tmp) / sum(exp(tmp)))
}


##########
library(DMwR)
for (i in 1:length(AthLeafAllPWMsImporatntPred)){
  AthLeafAllPWMsImporatntPred[[i]]$weight2 <- SoftMax(AthLeafAllPWMsImporatntPred[[i]]$weight)
}

l2 <- do.call(rbind,AthLeafAllPWMsImporatntPred)
total1 <- merge(df1,l2,by="PWMs")

ggplot(total1,aes(weight2,..density..))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))

ggplot(total1,aes(weight2,..density.., color=MotifCatalog))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))

################### zscore
for (i in 1:length(AthLeafAllPWMsImporatntPred)){
  AthLeafAllPWMsImporatntPred[[i]]$weight3 <- scale(AthLeafAllPWMsImporatntPred[[i]]$weight)
}

l2 <- do.call(rbind,AthLeafAllPWMsImporatntPred)
total1 <- merge(df1,l2,by="PWMs")

ggplot(total1,aes(weight3,..density..))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))

ggplot(total1,aes(weight3,..density.., color=MotifCatalog))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))





################### softmax2
for (i in 1:length(AthLeafAllPWMsImporatntPred)){
  AthLeafAllPWMsImporatntPred[[i]]$weight4 <- softmax2(AthLeafAllPWMsImporatntPred[[i]]$weight)
}

l2 <- do.call(rbind,AthLeafAllPWMsImporatntPred)
total1 <- merge(df1,l2,by="PWMs")

ggplot(total1,aes(weight4,..density..))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))

ggplot(total1,aes(weight4,..density.., color=MotifCatalog))  +
  geom_line(stat="density") + 
  theme_bw() + facet_wrap(.~Mark) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=10))



#####################
getcum <- function(total, markName){
  df <- subset(total,Mark== markName)
  df$count <- 1
  df1 <- df[order(df$weight,decreasing=T),] 
  l1 <- subset(df1, MotifCatalog == "Neutral motifs")
  l1$cumsum <- cumsum(l1$count)
  l1$cumsumPre <- l1$cumsum / sum(l1$count)
  
  l2 <- subset(df1, MotifCatalog == "GC-rich motif")
  l2$cumsum <- cumsum(l2$count)
  l2$cumsumPre <- l2$cumsum / sum(l2$count)
  
  l3 <- subset(df1, MotifCatalog == "AT-rich motif")
  l3$cumsum <- cumsum(l3$count)
  l3$cumsumPre <- l3$cumsum / sum(l3$count)
  
  l4 <- rbind(l1,l2,l3)
  l5 <- l4[order(l4$weight,decreasing=T),]
  l5$orders <- 1:nrow(l5)
  return(l5)
}

markNames <- names(AthLeafAllPWMsImporatntPred)
AthLeafAllPWMs_GC_cum <- lapply(1:length(AthLeafAllPWMsImporatntPred),
                                function(x) getcum(total = total,markName = markNames[x]))

names(AthLeafAllPWMs_GC_cum) <- markNames

AthLeafAllPWMs_GC_cum.df <- do.call(rbind,AthLeafAllPWMs_GC_cum)


AthLeafPWMsImportance_GC_countsPre <- ggplot(AthLeafAllPWMs_GC_cum.df,aes(x = orders,y=cumsumPre,color=MotifCatalog))+
  stat_ecdf()+facet_wrap(.~Mark) +
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# save : /home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafPWMsImportance_GC_countsPre.pdf





###################### 2021 7 7
getcum1 <- function(total, markName){
  Df1 <- subset(total,Mark== markName)
  Df1$count <- 1
  Df1 <- Df1[order(Df1$weight,decreasing=T),] 
  Df1$orders <- 1:nrow(Df1)
  Df2 <- Df1
  
  l1 <- Df2
  l1$count[which(l1$MotifCatalog != "Neutral motifs")] <- 0
  l1$cumsum <- cumsum(l1$count)
  Df1$Neutral <- 0
  Df1$Neutral[match(l1$PWMs,Df1$PWMs)] <- l1$cumsum
  Df1$NeutralRate <- Df1$Neutral / Df1$orders
  
  
  l2 <- Df2
  l2$count[which(l2$MotifCatalog != "GC-rich motif")] <- 0
  l2$cumsum <- cumsum(l2$count)
  Df1$GC <- 0
  Df1$GC[match(l2$PWMs,Df1$PWMs)] <- l2$cumsum
  Df1$GCRate <- Df1$GC / Df1$orders
  
  l3 <- Df2
  l3$count[which(l3$MotifCatalog != "AT-rich motif")] <- 0
  l3$cumsum <- cumsum(l3$count)
  Df1$AT <- 0
  Df1$AT[match(l3$PWMs,Df1$PWMs)] <- l3$cumsum
  Df1$ATRate <- Df1$AT / Df1$orders
  
  Df3 <- Df1[, c(1,8,10,12)]
  Df4 <- melt(Df3,id.vars = "PWMs")
  Df4$order <- rep(Df1$orders, times=3)
  Df4$Mark <- rep(Df1$Mark,times=3)
  return(Df4)
}

markNames <- names(AthLeafAllPWMsImporatntPred)
AthLeafAllPWMs_GC_cum <- lapply(1:length(AthLeafAllPWMsImporatntPred),
                                function(x) getcum1(total = total,markName = markNames[x]))

names(AthLeafAllPWMs_GC_cum) <- markNames

AthLeafAllPWMs_GC_cum.df <- do.call(rbind,AthLeafAllPWMs_GC_cum)


ggplot(data = AthLeafAllPWMs_GC_cum.df,aes(x=order,y=value,group = variable,color=variable))+
    geom_point(size=0.8)+
  geom_line()+facet_wrap(.~Mark) +
  xlab("PWMs")+#横坐标名称
  ylab("Rate")+#纵坐标名称
  theme_bw() 
  

## save : /home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafPWMsImportance_GC_countsPre_0707.pdf

# DO not need H2A
AthLeafAllPWMs_GC_cum_sub.df <- subset(AthLeafAllPWMs_GC_cum.df,Mark != "H2A")
ggplot(data = AthLeafAllPWMs_GC_cum_sub.df,aes(x=order,y=value,group = variable,color=variable))+
  geom_point(size=0.8)+
  geom_line()+facet_wrap(.~Mark) +
  xlab("PWMs")+#横坐标名称
  ylab("Rate")+#纵坐标名称
  theme_bw()


## save : /home/yanglab/Allwork/rework1/Ath/leaf/AthLeafPWM_reGC_AT/AthLeafPWMsImportance_GC_countsPre_0707_sub.pdf