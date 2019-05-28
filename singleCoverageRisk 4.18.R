# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)

#install.packages("sjPlot")
library(sjPlot)
library(grDevices)

#change color of CI fill here
fillcolor <- "#230FCB"
transparency <- .6

table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% SI','RR','95% SI')
table2$Effect <- c('Direct','Disseminated','Composite','Overall')
table4 <- table2

riskDiffRatio <- function(riskDF,tbl1,tbl2,coverage){
  
  naming2 <- paste(coverage,"Prev", sep = "") # name for prevalence
  naming4 <- paste(coverage,"Inc", sep = "") # name for incidence
  
# -----
  # individual prevalence

  individualDiff <- riskDF$risk11_prev - riskDF$risk01_prev
  tbl1[1,2] <-mean(individualDiff) #add to table

  #lowerPrevCI
  PrevCI<-quantile(individualDiff, probs = c(.025, .975))
  upperPrevCI <- PrevCI[2]
  lowerPrevCI <- PrevCI[1]
  tbl1[1,3] <- paste(lowerPrevCI,upperPrevCI) #find simulation intervals (2.5 and 97.5 percentiles)
  assign(paste("individualDiff",naming2,sep=""),c(mean(individualDiff),lowerPrevCI,upperPrevCI), envir = .GlobalEnv)
  individualRatio <- riskDF$risk11_prev/riskDF$risk01_prev
  tbl1[1,4] <- mean(individualRatio)
  CI<-quantile(individualRatio, probs = c(.025, .975))
  lowerCI <- CI[1]
  upperCI <- CI[2]
  tbl1[1,5] <- paste(lowerCI,upperCI)
  



  # incidence individual
  indIncDiff <-riskDF$risk11_inc - riskDF$risk01_inc
  tbl2[1,2] <- mean(indIncDiff) # risk difference for incidence
  #confidence intervals
  lowerIncCI <- quantile(indIncDiff, probs = .025)
  upperIncCI <- quantile(indIncDiff, probs = .975)
  tbl2[1,3] <- paste(lowerIncCI,upperIncCI)
  #assign(paste("individualDiff",naming4,sep=""),indIncDiff)
  assign(paste("fig4", coverage, sep = "_"), c(mean(indIncDiff),lowerIncCI,upperIncCI), envir = .GlobalEnv)
  indIncRatio <- riskDF$risk11_inc/riskDF$risk01_inc
  lowerIncCI <- quantile(indIncRatio, probs = .025)
  upperIncCI <- quantile(indIncRatio, probs = .975)
  tbl2[1,4] <- mean(indIncRatio)
  tbl2[1,5] <- paste(lowerIncCI,upperIncCI)

  #assign(paste("individualRatio",naming4,sep=""),c(mean(indIncRatio), lowerIncCI,upperIncCI), envir = .GlobalEnv)
  
  
  #-----
  #fix difference and ratios beyond this point
  # disseminated prevalence risk
  dissDiff <- riskDF$risk01_prev - riskDF$risk00_prev

  tbl1[2,2] <- mean(dissDiff)
  lowerCI <- quantile(dissDiff, probs = .025)
  upperCI <- quantile(dissDiff, probs = .975)
  tbl1[2,3] <- paste(lowerCI,upperCI)
  
  dissRatio <- riskDF$risk01_prev / riskDF$risk00_prev
  tbl1[2,4] <- mean(dissRatio)
  lowerCI <- quantile(dissRatio, probs = .025)
  upperCI <- quantile(dissRatio, probs = .975)
  tbl1[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming2, sep=""),dissDiff)
  assign(paste("disseminatedRatio",naming2, sep=""),dissRatio)
  
  #disseminated incident risk
  dissIncDiff <- riskDF$risk01_inc - riskDF$risk00_inc
  tbl2[2,2] <- mean(dissIncDiff)
  lowerCI <- quantile(dissIncDiff, probs = .025)
  upperCI <- quantile(dissIncDiff, probs = .975)
  tbl2[2,3] <- paste(lowerCI,upperCI)

  assign(paste("fig3", coverage, sep = "_"), c(mean(dissIncDiff),lowerCI,upperCI), envir = .GlobalEnv)
  
  dissIncRatio <- riskDF$risk01_inc / riskDF$risk00_inc
  tbl2[2,4] <- mean(dissIncRatio)
  lowerCI <- quantile(dissIncRatio, probs = .025)
  upperCI <- quantile(dissIncRatio, probs = .975)
  tbl2[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming4, sep=""),dissIncDiff)
  assign(paste("disseminatedRatio",naming4, sep=""),dissIncRatio)

  # overall prevalence risk
  
  overallDiff <- riskDF$risk1_prev - riskDF$risk00_prev
  tbl1[4,2] <- mean(overallDiff)
  lowerCI <- quantile(overallDiff, probs = .025)
  upperCI <- quantile(overallDiff, probs = .975)
  tbl1[4,3] <- paste(lowerCI,upperCI)
  
  overallRatio <-  riskDF$risk1_prev / riskDF$risk00_prev
  tbl1[4,4] <- mean(overallRatio)
  lowerCI <- quantile(overallRatio, probs = .025)
  upperCI <- quantile(overallRatio, probs = .975)
  tbl1[4,5] <- paste(lowerCI,upperCI)
  
  assign(paste("overallDiff",naming2, sep=""),overallDiff)
  assign(paste("overallRatio",naming2, sep=""),overallRatio)
  
  #overall incidence risk
  overallIncDiff <-  riskDF$risk1_inc - riskDF$risk00_inc
  tbl2[4,2] <- mean(overallIncDiff)
  lowerCI <- quantile(overallIncDiff, probs = .025)
  upperCI <- quantile(overallIncDiff, probs = .975)
  tbl2[4,3] <- paste(lowerCI,upperCI)
  overallIncRatio <- riskDF$risk1_inc / riskDF$risk00_inc
  tbl2[4,4] <- mean(overallIncRatio)
  lowerCI <- quantile(overallIncRatio, probs = .025)
  upperCI <- quantile(overallIncRatio, probs = .975)
  tbl2[4,5] <- paste(lowerCI,upperCI)
  
  assign(paste("overallDiff",naming4, sep=""),overallIncDiff)
  assign(paste("overallRatio",naming4, sep=""),overallIncRatio)
  
  # composite prevalence risk
  compositeDiff <- riskDF$risk11_prev - riskDF$risk00_prev
  tbl1[3,2] <- mean(compositeDiff)
  lowerCI <- quantile(compositeDiff, probs = .025)
  upperCI <- quantile(compositeDiff, probs = .975)
  tbl1[3,3] <- paste(lowerCI,upperCI)
  
  compositeRatio <- riskDF$risk11_prev / riskDF$risk00_prev
  tbl1[3,4] <- mean(compositeRatio)
  lowerCI <- quantile(compositeRatio, probs = .025)
  upperCI <- quantile(compositeRatio, probs = .975)
  tbl1[3,5] <- paste(lowerCI, upperCI)
  
  assign(paste("compositeDiff",naming2, sep=""),compositeDiff)
  assign(paste("overallRatio",naming2, sep=""),compositeRatio)
  
  # composite incidence risk
  compositeIncDiff <- riskDF$risk11_inc - riskDF$risk00_inc
  tbl2[3,2] <- mean(compositeIncDiff)
  lowerCI <- quantile(compositeIncDiff, probs = .025)
  upperCI <- quantile(compositeIncDiff, probs = .975)
  tbl2[3,3] <- paste(lowerCI,upperCI)
  
  compositeIncRatio <- riskDF$risk11_inc / riskDF$risk00_inc
  tbl2[3,4] <- mean(compositeIncRatio)
  lowerCI <- quantile(compositeIncRatio, probs = .025)
  upperCI <- quantile(compositeIncRatio, probs = .975)
  tbl2[3,5] <- paste(lowerCI,upperCI)
  assign(paste("compositeDiff",naming4, sep=""),compositeIncDiff)
  assign(paste("overallRatio",naming4, sep=""),compositeIncRatio)
  assign(paste("table2",coverage, sep="_"),tbl1, envir = .GlobalEnv)
  assign(paste("table4",coverage, sep="_"),tbl2, envir = .GlobalEnv)
}
riskDiffRatio(riskDF_10,table2, table4,10)
riskDiffRatio(riskDF_20,table2, table4,20)
riskDiffRatio(riskDF_30,table2, table4,30)
riskDiffRatio(riskDF_40,table2, table4,40)
riskDiffRatio(riskDF_50,table2, table4,50)
riskDiffRatio(riskDF_60,table2, table4,60)
riskDiffRatio(riskDF_70,table2, table4,70) 
riskDiffRatio(riskDF_80,table2, table4,80)
riskDiffRatio(riskDF_90,table2, table4,90)

for (i in 1:(length(listofcoverage))){riskDiffRatio(eval(as.name(listofcoverage[i])),table2,table4,listofcoverage[i])}
for (i in 1:length(listoftables)){
  thetable<-eval(as.name(listoftables[i]))
  flname<-paste(listoftables[i],".doc",sep="")
  x<-tab_df(x=thetable, show.rownames = F, alternate.rows = T, file = "x.doc")
  tab_df(eval(as.name(listoftables[1])),alternate.rows = T, file =paste(listoftables[i],".doc",sep=""))}
tab_df(eval(as.name(listoftables[31])),alternate.rows = T, file =paste(listoftables[31],".doc",sep=""))
#tables
setwd("/Volumes/GoogleDrive/My Drive/Networks/Brandon/Model Results/Tables")
#round table to 2 decimal places

tab_df(table2_30, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_30.doc")
tab_df(table4_30, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_30.doc")

tab_df(table2_70, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_70.doc")
tab_df(table4_70, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_70.doc")

tab_df(table4_40, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_40.doc")

tab_df(table2_10, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_10.doc")
tab_df(table4_10, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_10.doc")

tab_df(table2_20, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_20.doc")
tab_df(table4_20, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_20.doc")

tab_df(table2_40, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_40.doc")
tab_df(table4_40, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_40.doc")

tab_df(table2_50, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_50.doc")
tab_df(table4_50, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_50.doc")

tab_df(table2_60, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_60.doc")
tab_df(table4_60, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_60.doc")

tab_df(table2_80, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_80.doc")
tab_df(table4_80, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_80.doc")

tab_df(table2_90, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_90.doc")
tab_df(table4_90, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_90.doc")


#sensitivities
riskDiffRatio(riskactsDouble,table2,table4,"actsDouble")
riskDiffRatio(riskactsHalf,table2,table4,"actsHalf")
riskDiffRatio(riskrHalf,table2,table4,"rHalf")
riskDiffRatio(riskrDouble,table2,table4,"rDouble")

tab_df(table2_actsDouble, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_actsDouble.doc")
tab_df(table4_actsDouble, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_actsDouble.doc")

tab_df(table2_actsHalf, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_actsHalf.doc")
tab_df(table4_actsHalf, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_actsHalf.doc")

tab_df(table2_rHalf, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_rHalf.doc")
tab_df(table4_rHalf, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_rHalf.doc")

tab_df(table2_rDouble, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_rDouble.doc")
tab_df(table4_rDouble, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_rDouble.doc")

for (i in listofcoverage){
  for (j in listofcoverage){
    if (identical(eval(as.name(i)),eval(as.name(j)))){
      print(paste(i,j))
    }
  }
}
