# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)

install.packages("sjPlot")
library(sjPlot)
library(grDevices)

#change color of CI fill here
fillcolor <- "#230FCB"

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
  tbl1[1,2] <-round(mean(individualDiff),2) #add to table

  #lowerPrevCI
  PrevCI<-round(quantile(individualDiff, probs = c(.025, .975)),2)
  upperPrevCI <- PrevCI[2]
  lowerPrevCI <- PrevCI[1]
  tbl1[1,3] <- paste(lowerPrevCI,upperPrevCI) #find simulation intervals (2.5 and 97.5 percentiles)
  assign(paste("individualDiff",naming2,sep=""),c(mean(individualDiff),lowerPrevCI,upperPrevCI), envir = .GlobalEnv)
  individualRatio <- riskDF$risk11_prev/riskDF$risk01_prev
  tbl1[1,4] <- round(mean(individualRatio),2)
  CI<-round(quantile(individualRatio, probs = c(.025, .975)),2)
  lowerCI <- CI[1]
  upperCI <- CI[2]
  tbl1[1,5] <- paste(lowerCI,upperCI)
  



  # incidence individual
  indIncDiff <-riskDF$risk11_inc - riskDF$risk01_inc
  tbl2[1,2] <- round(mean(indIncDiff),2) # risk difference for incidence
  #confidence intervals
  lowerIncCI <- round(quantile(indIncDiff, probs = .025),2)
  upperIncCI <- round(quantile(indIncDiff, probs = .975),2)
  tbl2[1,3] <- paste(lowerIncCI,upperIncCI)
  #assign(paste("individualDiff",naming4,sep=""),indIncDiff)
  assign(paste("fig4", coverage, sep = "_"), c(mean(indIncDiff),lowerIncCI,upperIncCI), envir = .GlobalEnv)
  indIncRatio <- riskDF$risk11_inc/riskDF$risk01_inc
  lowerIncCI <- round(quantile(indIncRatio, probs = .025),2)
  upperIncCI <- round(quantile(indIncRatio, probs = .975),2)
  tbl2[1,4] <- round(mean(indIncRatio),2)
  tbl2[1,5] <- paste(lowerIncCI,upperIncCI)

  #assign(paste("individualRatio",naming4,sep=""),c(mean(indIncRatio), lowerIncCI,upperIncCI), envir = .GlobalEnv)
  
  
  #-----
  #fix difference and ratios beyond this point
  # disseminated prevalence risk
  dissDiff <- riskDF$risk01_prev - riskDF$risk00_prev

  tbl1[2,2] <- round(mean(dissDiff),2)
  lowerCI <- round(quantile(dissDiff, probs = .025),2)
  upperCI <- round(quantile(dissDiff, probs = .975),2)
  tbl1[2,3] <- paste(lowerCI,upperCI)
  
  dissRatio <- riskDF$risk01_prev / riskDF$risk00_prev
  tbl1[2,4] <- round(mean(dissRatio),2)
  lowerCI <- round(quantile(dissRatio, probs = .025),2)
  upperCI <- round(quantile(dissRatio, probs = .975),2)
  tbl1[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming2, sep=""),dissDiff)
  assign(paste("disseminatedRatio",naming2, sep=""),dissRatio)
  
  #disseminated incident risk
  dissIncDiff <- riskDF$risk01_inc - riskDF$risk00_inc
  tbl2[2,2] <- round(mean(dissIncDiff),2)
  lowerCI <- round(quantile(dissIncDiff, probs = .025),2)
  upperCI <- round(quantile(dissIncDiff, probs = .975),2)
  tbl2[2,3] <- paste(lowerCI,upperCI)

  assign(paste("fig3", coverage, sep = "_"), c(mean(dissIncDiff),lowerCI,upperCI), envir = .GlobalEnv)
  
  dissIncRatio <- riskDF$risk01_inc / riskDF$risk00_inc
  tbl2[2,4] <- round(mean(dissIncRatio),2)
  lowerCI <- round(quantile(dissIncRatio, probs = .025),2)
  upperCI <- round(quantile(dissIncRatio, probs = .975),2)
  tbl2[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming4, sep=""),dissIncDiff)
  assign(paste("disseminatedRatio",naming4, sep=""),dissIncRatio)

  # overall prevalence risk
  
  overallDiff <- riskDF$risk1_prev - riskDF$risk00_prev
  tbl1[4,2] <- round(mean(overallDiff),2)
  lowerCI <- round(quantile(overallDiff, probs = .025),2)
  upperCI <- round(quantile(overallDiff, probs = .975),2)
  tbl1[4,3] <- paste(lowerCI,upperCI)
  
  overallRatio <-  riskDF$risk1_prev / riskDF$risk00_prev
  tbl1[4,4] <- round(mean(overallRatio),2)
  lowerCI <- round(quantile(overallRatio, probs = .025),2)
  upperCI <- round(quantile(overallRatio, probs = .975),2)
  tbl1[4,5] <- paste(lowerCI,upperCI)
  
  assign(paste("overallDiff",naming2, sep=""),overallDiff)
  assign(paste("overallRatio",naming2, sep=""),overallRatio)
  
  #overall incidence risk
  overallIncDiff <-  riskDF$risk1_inc - riskDF$risk00_inc
  tbl2[4,2] <- round(mean(overallIncDiff),2)
  lowerCI <- round(quantile(overallIncDiff, probs = .025),2)
  upperCI <- round(quantile(overallIncDiff, probs = .975),2)
  tbl2[4,3] <- paste(lowerCI,upperCI)
  overallIncRatio <- riskDF$risk1_inc / riskDF$risk00_inc
  tbl2[4,4] <- round(mean(overallIncRatio),2)
  lowerCI <- round(quantile(overallIncRatio, probs = .025),2)
  upperCI <- round(quantile(overallIncRatio, probs = .975),2)
  tbl2[4,5] <- paste(lowerCI,upperCI)
  
  assign(paste("overallDiff",naming4, sep=""),overallIncDiff)
  assign(paste("overallRatio",naming4, sep=""),overallIncRatio)
  
  # composite prevalence risk
  compositeDiff <- riskDF$risk11_prev - riskDF$risk00_prev
  tbl1[3,2] <- round(mean(compositeDiff),2)
  lowerCI <- round(quantile(compositeDiff, probs = .025),2)
  upperCI <- round(quantile(compositeDiff, probs = .975),2)
  tbl1[3,3] <- paste(lowerCI,upperCI)
  
  compositeRatio <- riskDF$risk11_prev / riskDF$risk00_prev
  tbl1[3,4] <- round(mean(compositeRatio),2)
  lowerCI <- round(quantile(compositeRatio, probs = .025),2)
  upperCI <- round(quantile(compositeRatio, probs = .975),2)
  tbl1[3,5] <- paste(lowerCI, upperCI)
  
  assign(paste("compositeDiff",naming2, sep=""),compositeDiff)
  assign(paste("overallRatio",naming2, sep=""),compositeRatio)
  
  # composite incidence risk
  compositeIncDiff <- riskDF$risk11_inc - riskDF$risk00_inc
  tbl2[3,2] <- round(mean(compositeIncDiff),2)
  lowerCI <- round(quantile(compositeIncDiff, probs = .025),2)
  upperCI <- round(quantile(compositeIncDiff, probs = .975),2)
  tbl2[3,3] <- paste(lowerCI,upperCI)
  
  compositeIncRatio <- riskDF$risk11_inc / riskDF$risk00_inc
  tbl2[3,4] <- round(mean(compositeIncRatio),2)
  lowerCI <- round(quantile(compositeIncRatio, probs = .025),2)
  upperCI <- round(quantile(compositeIncRatio, probs = .975),2)
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



#tables
setwd("/Volumes/GoogleDrive/My Drive/Networks/Brandon/Model Results/Tables")
#round table to 2 decimal places

tab_df(table2_30, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_30.doc")
tab_df(table4_30, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_30.doc")

tab_df(table2_70, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_70.doc")
tab_df(table4_70, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_70.doc")



# figure 3 (incidence)
level <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
fig3Data <- as.data.frame(rbind(fig3_10,fig3_20,fig3_30,fig3_40,fig3_50,fig3_60,fig3_70,fig3_80, fig3_90))
fig3Data <- cbind(level, fig3Data)
colnames(fig3Data)<-c("level","estimate","lowerCI","upperCI")

fig4Data <- as.data.frame(rbind(fig4_10,fig4_20,fig4_30,fig4_40,fig4_50,fig4_60,fig4_70,fig4_80, fig4_90))
fig4Data <- cbind(level, fig4Data)
colnames(fig4Data)<-c("level","estimate","lowerCI","upperCI")

Iname <- expression(widehat(IE))
alph <- expression(alpha)
#plot figures 3 and 4
fig3plot <- ggplot(fig3Data) + geom_ribbon(aes(ymin=lowerCI, ymax=upperCI, x=level),alpha = 0.3, fill = fillcolor)+
  geom_line(aes(y=estimate, x=level))+
  xlab(expression(paste(alpha,"'",sep = "")))+
   ylab(expression(widehat(DE)(alpha)))

ggsave("Figure_3.eps",fig3plot)

fig4plot <- ggplot(fig4Data, aes(x = level, y= estimate)) + geom_line() +
  geom_ribbon(aes(ymin = lowerCI,ymax = upperCI, x = level), alpha = 0.3) +
  xlab(expression(paste(alpha,"'",sep = "")))+
  ylab(expression(widehat(IE)(alpha)))
ggsave("Figure_4.eps",fig4plot)
