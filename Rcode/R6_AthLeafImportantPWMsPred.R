######################################################################
# using different numbers of PWMs to predict histone modification levels 
# Redundant PWMs
# Data 2020 11 24
######################################################################################
AthLeafAllPWMPeakCountMatrix <-readRDS("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/results/AthLeafAllPWMPeakCountMatrix.RDS")


cv <- function(peakKmerRes.df, k=10)
{
  require(ranger)
  require(caret)
  folds <- createFolds(peakKmerRes.df$signal, k = k)
  lapply(folds, function(x)
  {
    cv_train <- peakKmerRes.df[-x, ]
    cv_test <- peakKmerRes.df[x, ]
    cv_model <- ranger(signal~., cv_train,num.threads = 4)
    cv_pred <- predict(cv_model, cv_test)
    corRes <- cor.test(cv_pred$prediction,cv_test$signal)
    return(c(corRes$p.value, corRes$estimate, cv_model$r.squared))
  }
  )
}

getPWMsSubPerformance <- function(markModel,j)
{
  library(dplyr)
  importantFeaturesOrder <- names(sort(markModel$variable.importance,decreasing = T))
  importantFeaturesOrderSub <- importantFeaturesOrder[1:j]
  PWMsSub_df <- select(PWMs_df,importantFeaturesOrderSub,"signal")
  PWMsSubPerformance <- cv(PWMsSub_df)
  PWMsSubSummary <- vector("list",length=2)
  names(PWMsSubSummary) <- c("markModelImp","PWMsSubPerformance")
  PWMsSubSummary$markModelImp <- importantFeaturesOrder
  PWMsSubSummary$PWMsSubPerformance <- PWMsSubPerformance
  return(PWMsSubSummary)
}


markNames <- names(AthLeafAllPWMPeakCountMatrix)


AthLeafAllPWMsImporatntPred <- vector("list",length = length(AthLeafAllPWMPeakCountMatrix))
names(AthLeafAllPWMsImporatntPred) <- names(AthLeafAllPWMPeakCountMatrix)


for (i in 1:length(AthLeafAllPWMPeakCountMatrix))
{
  aa <- as.numeric(floor(quantile(1:ncol(AthLeafAllPWMPeakCountMatrix[[i]])-1,probs=seq(0,1,0.05))))[-1]
  performances <- vector("list",length = length(aa))
  PWMs_df <- AthLeafAllPWMPeakCountMatrix[[i]]
  library(ranger)
  markModel <-ranger(signal~., PWMs_df,num.threads = 4,importance = "impurity")
  for (l in 1:length(aa))
  {
    performances[[l]] <- getPWMsSubPerformance(markModel=markModel,j=aa[l])
  }
  AthLeafAllPWMsImporatntPred[[i]] <- performances
  print(sprintf("i=%d",i))
}




saveRDS(AthLeafAllPWMsImporatntPred,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/AthLeafPWMsImorptancePred/AthLeafImportantPWMsPredPerformace.RDS")



getEachMarkperformance <- function(AthLeafEachPWMsImporatntPred,ll)
{
  aaa <- data.frame(matrix(ncol = 5,nrow = 40))
  colnames(aaa) <- c("Marks","PWMsNum","Type","Performance","sd")
  aaa$Marks <- rep(markNames[ll],40)
  aaa$PWMsNum <- rep(as.numeric(floor(quantile(1:ncol(AthLeafAllPWMPeakCountMatrix[[ll]])-1,probs=seq(0,1,0.05))))[-1],2)
  aaa$Type <- rep(c("PCC","Rsquared"),each=20)
  for (j in 1:20) 
  {
    j20 <- j + 20
    aaa$Performance[j] <- mean(unlist(lapply(AthLeafEachPWMsImporatntPred[[j]]$PWMsSubPerformance, function(x) as.numeric(x[2]))))
    aaa$Performance[j20] <- mean(unlist(lapply(AthLeafEachPWMsImporatntPred[[j]]$PWMsSubPerformance, function(x) as.numeric(x[3]))))
    aaa$sd[j] <- sd(unlist(lapply(AthLeafEachPWMsImporatntPred[[j]]$PWMsSubPerformance, function(x) as.numeric(x[2]))))
    aaa$sd[j20] <- sd(unlist(lapply(AthLeafEachPWMsImporatntPred[[j]]$PWMsSubPerformance, function(x) as.numeric(x[3]))))
    
  }
  return(aaa)
}




AthLeafAllPWMsImporatntPred_PCC_R <- vector("list",length = length(AthLeafAllPWMPeakCountMatrix))
names(AthLeafAllPWMsImporatntPred_PCC_R) <- names(AthLeafAllPWMPeakCountMatrix)

for (i in 1:length(AthLeafAllPWMsImporatntPred_PCC_R))
{
  AthLeafAllPWMsImporatntPred_PCC_R[[i]] <- getEachMarkperformance(AthLeafEachPWMsImporatntPred = AthLeafAllPWMsImporatntPred[[i]],
                                                                   ll = i)
}

AthLeafAllPWMsImporatntPred_PCC_R.df <- do.call("rbind",AthLeafAllPWMsImporatntPred_PCC_R)
saveRDS(AthLeafAllPWMsImporatntPred_PCC_R.df,
        file = "/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/AthLeafPWMsImorptancePred/AthLeafAllPWMsImporatntPred_PCC_R_df.RDS")
library(ggplot2)


ggplot(AthLeafAllPWMsImporatntPred_PCC_R.df, 
       aes(x=PWMsNum, y=Performance, color=Type,fill=Type)) + geom_line(size=0.5) + 
  facet_wrap(.~Marks, nrow = 3, scales = "free")+
  theme_bw() + # whithe
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0,130),n.breaks = 20,expand=c(0,0)) +
  geom_errorbar(aes(ymin=Performance-sd, ymax=Performance+sd), width=.2,color="black",
                position=position_dodge(1))




#####################################################################################################
####################################

######################################################################################################
# get the best information PWMs from top
# 2020 11 27
##################################################

AthLeafAllPWMs1_mis1 <- readRDS("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/NEWet1.0Mis1/AthLeafAllPWMs1_mis1.RDS")
AthLeafAllPWMsImporatntPred_PCC_R.df <- readRDS("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/AthLeafPWMsImorptancePred/AthLeafAllPWMsImporatntPred_PCC_R_df.RDS")
markNames <- names(AthLeafAllPWMs1_mis1)

getEachMarkImportantPWMs <- function(AthLeafAllPWMs1_mis1,AthAthLeafAllPWMsImporatntPred_PCC_R.df,i,AthLeafAllPWMsImporatntPred)
{
  markPerformance <- subset(AthLeafAllPWMsImporatntPred_PCC_R.df,AthLeafAllPWMsImporatntPred_PCC_R.df$Marks == markNames[i])
  bestPWMsNum <- markPerformance$PWMsNum[which(diff(markPerformance$Performance[1:20]) < 0.005)[1]]
  bestPWMsNames <- AthLeafAllPWMsImporatntPred[[i]][[1]]$markModelImp[1:bestPWMsNum]
  names(AthLeafAllPWMs1_mis1[[i]]) <- unlist(lapply(AthLeafAllPWMs1_mis1[[i]], function(x) x@name))
  bestPWMs <- AthLeafAllPWMs1_mis1[[i]][names(AthLeafAllPWMs1_mis1[[i]]) %in% bestPWMsNames]
  return(bestPWMs)
}


AthLeafAllPWMPeakCountMatrix <- readRDS("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/results/AthLeafAllPWMPeakCountMatrix.RDS")

AthLeafAllPWMsImporatntPred  <- readRDS("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/AthLeafPWMsImorptancePred/AthLeafImportantPWMsPredPerformace.RDS")

getZeroPWMs_df <- function(AthLeafAllPWMsImporatntPred,AthLeafAllPWMPeakCountMatrix,featureCut,i)
{
  library(dplyr)
  AllFeature_df <- AthLeafAllPWMPeakCountMatrix[[i]]
  Features <- AthLeafAllPWMsImporatntPred[[i]][[1]]$markModelImp
  FeatureTop <- Features[1:featureCut]
  FeatureTop_df <- AllFeature_df %>% select(FeatureTop)
  FeatureLast <-Features[(featureCut+1):length(Features)]
  FeatureLast_df <- AllFeature_df %>% select(FeatureLast)
  FeatureLast_df[FeatureLast_df >= 0] =0
  FeatureLast_Top_df <- cbind(FeatureTop_df,FeatureLast_df)
  FeatureLast_Top_df$signal <- AllFeature_df$signal
  FeatureLast_Top_results <- cv(FeatureLast_Top_df)
}

for (i in 1:length(h3k4me3PWMs))
{
  names(h3k4me3PWMs)[i] <- h3k4me3PWMs[[i]]@name
}



############## 2021 7 8
# exclude H2A
AthLeafAllPWMsImporatntPred_PCC_R.df <- readRDS("/home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/AthLeafPWMsImorptancePred/AthLeafAllPWMsImporatntPred_PCC_R_df.RDS")
library(ggplot2)
AthLeafAllPWMsImporatntPred_PCC_R_sub.df <- subset(AthLeafAllPWMsImporatntPred_PCC_R.df,Marks != "H2A")
Rsquared_95 <-c(24,23,25,15,11,23,12,13,22,25,25,28)

AthLeafAllPWMsImporatntPred_PCC_R_sub.df$Rsquared_95 <- rep(Rsquared_95,each=40)
ggplot(AthLeafAllPWMsImporatntPred_PCC_R_sub.df, 
       aes(x=PWMsNum, y=Performance, color=Type,fill=Type)) + geom_line(size=0.5) +
  geom_vline(aes(xintercept = Rsquared_95), color="red",
             linetype="dashed", size=1) + 
  facet_wrap(.~Marks, nrow = 3, scales = "free")+
  theme_bw() + # whithe
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0,130),n.breaks = 20,expand=c(0,0)) +
  geom_errorbar(aes(ymin=Performance-sd, ymax=Performance+sd), width=.2,color="black",
                position=position_dodge(1))
# save : /home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/AthLeafPWMsImorptancePred/AthLeafAllPWMsImporatntPred_PCC_R_df_sub.pdf


############# 2021 8 18 at home
h3k4me3 <- subset(AthLeafAllPWMsImporatntPred_PCC_R_sub.df,Marks=="H3K4me3")
ggplot(h3k4me3, 
       aes(x=PWMsNum, y=Performance, color=Type,fill=Type)) + geom_line(size=0.7) +
  geom_vline(aes(xintercept = Rsquared_95), color="red",
             linetype="dashed", size=1) + 
  theme_bw() + # whithe
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0,130),n.breaks = 20,expand=c(0,0)) +
  geom_errorbar(aes(ymin=Performance-sd, ymax=Performance+sd), width=.5,color="black",
                position=position_dodge(1))
# save /home/yanglab/Allwork/reWork/part2/AthLeafPWMsPredMis1/AthLeafPWMsImorptancePred/AthLeaf_h3k4me3_ImporatntPred_PCC_R_df_sub
# 4 * 5









