# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)
## notes: need to make separate function for all estimators (include quantiles)
## then: through this, do risk calculations
# also need to add all to separate table for plotting?
library(sjPlot)
#have to have cum and prev ave first
table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% CI','RR','95% CI')
table2$Effect <- c('Individual','Disseminated','Composite','Overall')
table4 <- table2

riskDiffRatio <- function(riskDF,tbl1,tbl2,coverage){ #give the file of endtime data
  # and prevalence or incidence
  # may move the risks down when using estimators
  
  naming2 <- paste(coverage,"Prev", sep = "") # name for prevalence
  naming4 <- paste(coverage,"Inc", sep = "") # name for incidence
  
# -----
  # individual prevalence
  #individualDiff <- mean(riskDF$risk11) - mean(riskDF$risk01)
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
  
  #assign(paste("individualRatio",naming2,sep=""),mean(individualRatio), lowerCI, upperCI)


  # incidence individual
  indIncDiff <-riskDF$risk11_inc - riskDF$risk01_inc
  tbl2[1,2] <- mean(indIncDiff) # risk difference for incidence
  #confidence intervals
  lowerIncCI <- quantile(indIncDiff, probs = .025)
  upperIncCI <- quantile(indIncDiff, probs = .975)
  tbl2[1,3] <- paste(lowerIncCI,upperIncCI)
  #assign(paste("individualDiff",naming4,sep=""),indIncDiff)
  assign(paste("fig3", coverage, sep = "_"), c(mean(indIncDiff),lowerIncCI,upperIncCI), envir = .GlobalEnv)
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
  assign(paste("fig4", coverage, sep = "_"), c(dissIncDiff,lowerCI,upperCI), envir = .GlobalEnv)
  
  dissIncRatio <- riskDF$risk01_inc / riskDF$risk00_inc
  tbl2[2,4] <- mean(dissIncRatio)
  lowerCI <- quantile(dissIncRatio, probs = .025)
  upperCI <- quantile(dissIncDiff, probs = .975)
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

riskDiffRatio(risk30,table2, table4,30) 
riskDiffRatio(risk70,table2, table4,70) 

#tables
tab_df(table2_30, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_30.doc")
tab_df(table4_30, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_30.doc")

tab_df(table2_70, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2_70.doc")
tab_df(table4_70, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4_70.doc")



# figure 3 (incidence)
level <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
fig3Data <- as.data.frame(rbind(fig3_10,fig3_20,fig3_30,fig3_40,fig3_50,fig3_60,fig3_70,fig3_80, fig3_90))
fig3Data <- cbind(level, fig3Data)
colnames(fig3Data)<-c("level","estimate","lowerCI","upperCI")


#plot figures 3 and 4
fig3plot <- ggplot(fig3Data) + geom_ribbon(aes(ymin=lowerCI, ymax=upperCI, x=level),alpha = 0.3)+geom_line(aes(y=estimate, x=level))+
  scale_colour_manual("",values="blue")+scale_fill_manual("",values="blue")

ggsave("Figure_3.eps",fig3plot)

fig4plot <- ggplot(fig4Data, aes(x = level, y= estimate)) + geom_line()
  geom_ribbon(aes(ymin = lowerCI,ymax = upperCI))
ggsave("Figure_4.eps",fig4plot)
