# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)
## notes: need to make separate function for all estimators (include quantiles)
## then: through this, do risk calculations
# also need to add all to separate table for plotting?

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
  
  #quantiles
  Prev_11_quants <- quantile(riskDF$risk11_prev,probs = c(.025,.975),na.rm = T)
  Prev_01_quants <- quantile(riskDF$risk01_prev,probs = c(.025,.975),na.rm = T)
  Prev_00_quants <- quantile(riskDF$risk00_prev,probs = c(.025,.975),na.rm = T)
  Prev_1_quants <- quantile(riskDF$risk1_prev,probs = c(.025,.975),na.rm = T)
  
  Inc_11_quants <- quantile(riskDF$risk11_inc,probs = c(.025,.975),na.rm = T)
  Inc_01_quants <- quantile(riskDF$risk01_inc,probs = c(.025,.975),na.rm = T)
  Inc_00_quants <- quantile(riskDF$risk00_inc,probs = c(.025,.975),na.rm = T)
  Inc_1_quants <- quantile(riskDF$risk1_inc,probs = c(.025,.975),na.rm = T)

  # create differences and ratios for both prevalence and incidence
  # individual prevalence
  #individualDiff <- mean(riskDF$risk11) - mean(riskDF$risk01)
  individualDiff <- mean(riskDF$risk11_prev) - mean(riskDF$risk01_prev)
  tbl1[1,2] <- individualDiff #add to table


  #lowerPrevCI
  lowerPrevCI<-Prev_11_quants[1]-Prev_01_quants[2]
  upperPrevCI <- Prev_11_quants[2]-Prev_01_quants[1]
  tbl1[1,3] <- paste(lowerPrevCI,upperPrevCI) #find simulation intervals (2.5 and 97.5 percentiles)
  assign(paste("individualDiff",naming2,sep=""),c(individualDiff,lowerPrevCI,upperPrevCI), envir = .GlobalEnv)
  individualRatio <- mean(riskDF$risk11_prev)/mean(riskDF$risk01_prev)
  tbl1[1,4] <- individualRatio
  lowerCI <- Prev_11_quants[1]/Prev_01_quants[2]
  upperCI <- Prev_11_quants[2]/Prev_01_quants[1]
  tbl1[1,5] <- paste(lowerCI,upperCI)
  
  assign(paste("individualRatio",naming2,sep=""),individualRatio)


  # incidence individual
  indIncDiff <- mean(riskDF$risk11_inc)-mean(riskDF$risk01_inc)
  tbl2[1,2] <- indIncDiff # risk difference for incidence
  #confidence intervals
  lowerIncCI <- Inc_11_quants[1]-Inc_01_quants[2]
  upperIncCI <- Inc_11_quants[2]-Inc_01_quants[1]
  tbl2[1,3] <- paste(lowerIncCI,upperIncCI)
  assign(paste("fig3", coverage, sep = "_"), c(indIncDiff,lowerIncCI,upperIncCI), envir = .GlobalEnv)
  indIncRatio <- mean(riskDF$risk11_inc)/mean(riskDF$risk01_inc)
  lowerIncCI <- Inc_11_quants[1]/Inc_01_quants[2]
  upperIncCI <- Inc_11_quants[2]/Inc_01_quants[1]
  tbl2[1,4] <- indIncRatio
  tbl2[1,5] <- paste(lowerIncCI,upperIncCI)
  assign(paste("individualDiff",naming4,sep=""),indIncDiff)
  assign(paste("individualRatio",naming4,sep=""),indIncRatio)
  
  
  #-----
  #fix difference and ratios beyond this point
  # disseminated prevalence risk
  dissDiff <- mean(riskDF$risk01_prev) - mean(riskDF$risk00_prev)

  tbl1[2,2] <- dissDiff
  lowerCI <- Prev_01_quants[1] - Prev_00_quants[2]
  upperCI <- Prev_01_quants[2] - Prev_00_quants[1]
  tbl1[2,3] <- paste(lowerCI,upperCI)
  
  dissRatio <- mean(riskDF$risk01_prev) / mean(riskDF$risk00_prev)
  tbl1[2,4] <- dissRatio
  lowerCI <- Prev_01_quants[1] / Prev_00_quants[2]
  upperCI <- Prev_01_quants[2] / Prev_00_quants[1]
  tbl1[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming2, sep=""),dissDiff)
  assign(paste("disseminatedRatio",naming2, sep=""),dissRatio)
  
  #disseminated incident risk
  dissIncDiff <- mean(riskDF$risk01_inc) - mean(riskDF$risk00_inc)
  tbl2[2,2] <- dissIncDiff
  lowerCI <- Inc_01_quants[1] - Inc_00_quants[2]
  upperCI <- Inc_01_quants[2] - Inc_00_quants[1]
  tbl2[2,3] <- paste(lowerCI,upperCI)
  assign(paste("fig4", coverage, sep = "_"), c(dissIncDiff,lowerCI,upperCI), envir = .GlobalEnv)
  
  dissIncRatio <- mean(riskDF$risk01_inc) / mean(riskDF$risk00_inc)
  tbl2[2,4] <-dissIncRatio
  lowerCI <- Inc_01_quants[1] / Inc_00_quants[2]
  upperCI <- Inc_01_quants[2] / Inc_00_quants[1]
  tbl2[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming4, sep=""),dissIncDiff)
  assign(paste("disseminatedRatio",naming4, sep=""),dissIncRatio)

  # overall prevalence risk
  
  overallDiff <- mean(riskDF$risk1_prev) - mean(riskDF$risk00_prev)
  tbl1[4,2] <- overallDiff
  lowerCI <- Prev_1_quants[1] - Prev_00_quants[2]
  upperCI <- Prev_1_quants[2] - Prev_00_quants[1]
  tbl1[4,3] <- paste(lowerCI,upperCI)
  
  overallRatio <-  mean(riskDF$risk1_prev) / mean(riskDF$risk00_prev)
  tbl1[4,4] <- overallRatio
  lowerCI <- Prev_1_quants[1] / Prev_00_quants[2]
  upperCI <- Prev_1_quants[2] / Prev_00_quants[1]
  tbl1[4,5] <- paste(lowerCI,upperCI)
  
  assign(paste("overallDiff",naming2, sep=""),overallDiff)
  assign(paste("overallRatio",naming2, sep=""),overallRatio)
  
  #overall incidence risk
  overallIncDiff <-  mean(riskDF$risk1_inc) - mean(riskDF$risk00_inc)
  tbl2[4,2] <- overallIncDiff
  lowerCI <- Inc_1_quants[1] -Inc_00_quants[2]
  upperCI <- Inc_1_quants[2] - Inc_00_quants[1]
  tbl2[4,3] <- paste(lowerCI,upperCI)
  overallIncRatio <- mean(riskDF$risk1_inc) / mean(riskDF$risk00_inc)
  tbl2[4,4] <- overallIncRatio
  lowerCI <- Inc_1_quants[1] / Inc_00_quants[2]
  upperCI <- Inc_1_quants[2] / Inc_00_quants[1]
  tbl2[4,5] <- paste(lowerCI,upperCI)
  
  assign(paste("overallDiff",naming4, sep=""),overallIncDiff)
  assign(paste("overallRatio",naming4, sep=""),overallIncRatio)
  
  # composite prevalence risk
  compositeDiff <- mean(riskDF$risk11_prev) - mean(riskDF$risk00_prev)
  tbl1[3,2] <- compositeDiff
  lowerCI <- Inc_11_quants[1] -Inc_00_quants[2]
  upperCI <- Inc_11_quants[2] - Inc_00_quants[1]
  tbl1[3,3] <- paste(lowerCI,upperCI)
  
  compositeRatio <- mean(riskDF$risk11_prev) / mean(riskDF$risk00_prev)
  tbl1[3,4] <- compositeRatio
  lowerCI <- Prev_11_quants[1] /Prev_00_quants[2]
  upperCI <- Prev_11_quants[2] / Prev_00_quants[1]
  tbl1[3,5] <- paste(lowerCI, upperCI)
  
  assign(paste("compositeDiff",naming2, sep=""),compositeDiff)
  assign(paste("overallRatio",naming2, sep=""),compositeRatio)
  
  # composite incidence risk
  compositeIncDiff <- mean(riskDF$risk11_inc) - mean(riskDF$risk00_inc)
  tbl2[3,2] <- compositeIncDiff
  lowerCI <- Inc_11_quants[1] -Inc_00_quants[2]
  upperCI <- Inc_11_quants[2] - Inc_00_quants[1]
  tbl2[3,3] <- paste(lowerCI,upperCI)
  compositeIncRatio <- mean(riskDF$risk11_inc) / mean(riskDF$risk00_inc)
  tbl2[3,4] <- compositeIncRatio
  lowerCI <- Inc_11_quants[1] /Inc_00_quants[2]
  upperCI <- Inc_11_quants[2] / Inc_00_quants[1]
  tbl2[3,5] <- paste(lowerCI,upperCI)
  assign(paste("compositeDiff",naming4, sep=""),compositeIncDiff)
  assign(paste("overallRatio",naming4, sep=""),compositeIncRatio)
  assign(paste("table2",coverage, sep="_"),tbl1, envir = .GlobalEnv)
  assign(paste("table4",coverage, sep="_"),tbl2, envir = .GlobalEnv)
}

riskDiffRatio(risk10,table2, table4,10) #treatment 10%


#tables
tab_df(table2_50, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2.doc")
tab_df(table4_50, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4.doc")


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
