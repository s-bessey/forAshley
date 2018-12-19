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
  #may move the risks down when using estimators
  naming2 <- paste(coverage,"Prev", sep = "") # name for prevalence
  naming4 <- paste(coverage,"Inc", sep = "") # name for incidence
  
  #quantiles
  Prev_11_quants <- quantile(riskDF$risk11_prev,probs = .025,.975)
  Prev_01_quants <- quantile(riskDF$risk01_prev,probs = .025,.975)
  Prev_00_quants <- quantile(riskDF$risk00_prev,probs = .025,.975)
  
  Inc_11_quants <- quantile(riskDF$risk11_inc,probs = .025,.975)
  Inc_01_quants <- quantile(riskDF$risk01_inc,probs = .025,.975)
  Inc_00_quants <- quantile(riskDF$risk00_inc,probs = .025,.975)
  
  # create differences and ratios for both prevalence and incidence
  # individual prevalence
  #individualDiff <- mean(riskDF$risk11) - mean(riskDF$risk01)
  individualDiff <- riskDF$risk11 - riskDF$risk01
  tbl1[1,2] <- individualDiff #add to table
  #lowerPrevCI
  lowerPrevCI<-Prev_11_quants[1]-Prev_01_quants[2]
  upperPrevCI <- Prev_11_quants[2]-Prev_01_quants[1]
  tbl1[1,3] <- paste(lowerPrevCI,upperPrevCI) #find simulation intervals (2.5 and 97.5 percentiles)
  
  individualRatio <- mean(riskDF$risk11_prev)/mean(riskDF$risk01_prev)
  tbl1[1,4] <- individualRatio
  lowerCI <- Prev_11_quants[1]/Prev_01_quants[2]
  upperCI <- Prev_11_quants[2]/Prev_01_quants[1]
  tbl1[1,5] <- paste(lowerCI,upperCI)
  assign(paste("individualDiff",naming2,sep=""),individualDiff)
  assign(paste("individualRatio",naming2,sep=""),individualRatio)


  # incidence individual
  indIncDiff <- mean(riskDF$risk11_inc)-mean(riskDF$risk01_inc)
  tbl2[1,2] <- indIncDiff # risk difference for incidence
  #confidence intervals
  lowerIncCI <- Inc_11_quants[1]-Inc_01_quants[2]
  upperIncCI <- Inc_11_quants[2]-Inc_01_quants[1]
  tbl2[1,3] <- paste(lowerIncCI,upperIncCI)
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

  tbl1[2,2] <- dissDiff$estimate
  lowerCI <- Prev_01_quants[1] - Prev_00_quants[2]
  upperCI <- Prev_01_quants[2] - Prev_00_quants[1]
  tbl1[2,3] <- paste(lowerCI,upperCI)
  
  dissRatio <- mean(riskDF$risk01_prev) / mean(riskDF$risk00_prev)
  tbl1[2,4] <- dissRatio$estimate
  lowerCI <- Prev_01_quants[1] / Prev_00_quants[2]
  upperCI <- Prev_01_quants[2] / Prev_00_quants[1]
  tbl1[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming2, sep=""),dissDiff)
  assign(paste("disseminatedRatio",naming2, sep=""),dissRatio)
  
  #disseminated incident risk
  dissIncDiff <- mean(riskDF$risk01_inc) - mean(riskDF$risk00_inc)
  tbl2[2,2] <- mean(dissIncDiff)[1]
  lowerCI <- Inc_01_quants[1] - Inc_00_quants[2]
  upperCI <- Inc_01_quants[2] - Inc_00_quants[1]
  tbl2[2,3] <- paste(lowerCI,upperCI)
  
  dissIncRatio <- mean(riskDF$risk01_inc) / mean(riskDF$risk00_inc)
  tbl2[2,4] <-dissIncRatio
  lowerCI <- Inc_01_quants[1] / Inc_00_quants[2]
  upperCI <- Inc_01_quants[2] / Inc_00_quants[1]
  tbl2[2,5] <- paste(lowerCI,upperCI)
  
  assign(paste("disseminatedDiff",naming4, sep=""),dissIncDiff)
  assign(paste("disseminatedRatio",naming4, sep=""),dissIncRatio)

  # overall prevalence risk
  
  overallDiff <- mean(riskDF$risk1_prev) - mean(riskDF$risk0_prev)
  tbl1[4,2] <- overallDiff
  lowerCI <- Prev_1_quants[1] - Prev_0_quants[2]
  upperCI <- Prev_1_quants[2] - Prev_0_quants[1]
  tbl1[4,3] <- paste(lowerCI,upperCI)
  
  overallRatio <-  mean(riskDF$risk1_prev) / mean(riskDF$risk0_prev)
  tbl1[4,4] <- overallRatio
  lowerCI <- Prev_1_quants[1] / Prev_0_quants[2]
  upperCI <- Prev_1_quants[2] / Prev_0_quants[1]
  tbl1[4,5] <- paste(lowerCI,upperCI)
  
  assign(print("overallDiff",naming2, sep=""),overallDiff)
  assign(print("overallRatio",naming2, sep=""),overallRatio)
  
  #overall incidence risk
  overallIncDiff <-  mean(riskDF$risk1_inc) - mean(riskDF$risk0_inc)
  tbl2[4,2] <- overallIncDiff
  lowerCI <- Inc_1_quants[1] -Inc_0_quants[2]
  upperCI <- Inc_1_quants[2] - Inc_0_quants[1]
  tbl2[4,3] <- paste(lowerCI,upperCI)
  overallIncRatio <- mean(riskDF$risk1_inc) / mean(riskDF$risk0_inc)
  tbl2[4,4] <- overallIncRatio
  lowerCI <- Inc_1_quants[1] / Inc_0_quants[2]
  upperCI <- Inc_1_quants[2] / Inc_0_quants[1]
  tbl2[4,5] <- paste(lowerCI,upperCI)
  
  assign(print("overallDiff",naming4, sep=""),overallIncDiff)
  assign(print("overallRatio",naming4, sep=""),overallIncRatio)
  
  # composite prevalence risk
  compositeDiff <- mean(riskDF$risk11_prev) - mean(riskDF$risk00_prev)
  tbl1[3,2] <- compositeDiff$estimate
  lowerCI <- Inc_11_quants[1] -Inc_00_quants[2]
  upperCI <- Inc_11_quants[2] - Inc_00_quants[1]
  tbl1[3,3] <- paste(lowerCI,upperCI)
  
  compositeRatio <- mean(riskDF$risk11_prev) / mean(riskDF$risk00_prev)
  tbl1[3,4] <- mean(compositeRatio$estimate)[1]
  tbl1[3,5] <- paste(quantile(compositeRatio$estimate, probs = .025),
                     quantile(compositeRatio$estimate, probs = .975))
  
  assign(print("compositeDiff",naming2, sep=""),compositeDiff)
  assign(print("overallRatio",naming2, sep=""),compositeRatio)
  
  # composite incidence risk
  compositeIncDiff <- mean(riskDF$risk11_inc) - mean(riskDF$risk00_inc)
  compositeRatio <- mean(riskDF$risk11_inc) / mean(riskDF$risk00_inc)
  assign(print("compositeDiff",naming4, sep=""),compositeIncDiff)
  assign(print("overallRatio",naming4, sep=""),compositeIncRatio)
}

riskDiffRatio() #treatment 10%


#tables
tab_df(table2, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2.]\doc")
tab_df(table4, col.header = colnames(table4), show.rownames = F, alternate.rows = T,file="Table4.]\doc")


# figure 3 (incidence)
fig3Data <- as.data.frame(matrix(nrow = 10, ncol = 4))
colnames(fig3Data) <- c("level", "estimate","lowerCI","upperCI")
fig4Data <- fig3Data

figure34dataCreation <- function(treatment){
  riskDiffRatio(treatment)
  level <- substr(deparse(substitute(treatment)), nchar(treatment)-1,nchar(treatment))
  figure3Data[level,1] <- paste('.', level,sep = "")
  figure3Data[level,2] <- indIncDiff
  figure3Data[level,3] <- lowerIncCI
  figure3Data[level,1] <- upperIncCI
  
  figure4Data[level,1] <- paste('.', level,sep = "")
  figure4Data[level,2] <- indPrevDiff
  figure4Data[level,3] <- lowerPrevCI
  figure4Data[level,4] <- upperPrevCI
}

#plot figures 3 and 4
fig3plot <- ggplot(fig3Data, aes(x = level, y= estimate)) +
  geom_ribbon(aes(ymin = lowerCI,ymax = upperCI))
ggsave("Figure_3.eps",fig3plot)

fig4plot <- ggplot(fig4Data, aes(x = level, y= estimate)) +
  geom_ribbon(aes(ymin = lowerCI,ymax = upperCI))
ggsave("Figure_4.eps",fig4plot)
