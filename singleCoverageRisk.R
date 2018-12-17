# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)


#have to have cum and prev ave first
table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% CI','RR','95% CI')
table2$Effect <- c('Individual','Disseminated','Composite','Overall')
table4 <- table2

riskDiffRatio <- function(treatment,tbl1,tbl2){ #give the file of endtime data
  # and prevalence or incidence
  
  naming2 <- paste(treatment$coverage[1],"Prev", sep = "") # name for prevalence
  naming4 <- paste(treatment$coverage[1],"Inc", sep = "") # name for incidence
  
  N_11 <- treatment$N_PrEP[treatment$t == max(treatment$t)] # number of people on PrEP
  HIV_11 <- treatment$N_HIV_PrEP[treatment$t == max(treatment$t)] # prevalence of HIV+ PrEP-users
  Inc_11 <- treatment$N_HIV_PrEP[treatment$t == max(treatment$t)] -
    treatment$N_HIV_PrEP[treatment$t == 0] #cum incidence of HIV in PrEP users
    
  N_01 <- treatment$N_noPrEP[treatment$t == max(treatment$t)] # number of people not on PrEP in treatment gp
  HIV_01 <- treatment$HIV_noPreP[treatment$t == max(treatment$t)] # HIV+ people not on PrEP in treatment gp
  Inc_01 <- treatment$Inc_HIV_noPrEP[treatment$t == max(treatment$t)] -
    treatment$N_HIV_PrEP[treatment$t == 0] # HIV incidence in no-PrEP treatment gp
  
  N_00 <- treatment$N_control[treatment$t == max(treatment$t)] # population of control gp
  HIV_00 <- treatment$HIV_control[treatment$t == max(treatment$t)] # HIV prevalence in control
  Inc_00 <- treatment$Inc_control[treatment$t == max(treatment$t)] -
    treatment$N_HIV_PrEP[treatment$t == 0] # HIV incidence in control
  
  
  # create differences and ratios for both prevalence and incidence
  # individual prevalence
  individualDiff <- riskdifference(HIV_11,HIV_01,# individual risk difference
                             N_11,N_01)
  tbl1[1,2] <- ave(individualDiff$estimate)[1] #add to table
  assign(paste("lowerPrevCI"), quantile(individualIncDiff$estimate, probs = .025))
  assign(paste("upperPrevCI"),quantile(individualIncDiff$estimate, probs = .975))
  tbl1[1,3] <- paste(quantile(individualDiff$estimate, probs = .025),
                     quantile(individualDiff$estimate, probs = .975)) #find simulation intervals (2.5 and 97.5 percentiles)
  individualRatio <- riskratio(HIV_11,HIV_01, # individual risk ratio
                               N_11,N_01)
  tbl1[1,4] <- ave(individualRatio$estimate)[1]
  tbl1[1,5] <- paste(quantile(individualRatio$estimate, probs = .025),
                     quantile(individualRatio$estimate, probs = .975))
  assign(paste("individualDiff",naming2,sep=""),individualDiff)
  assign(paste("individualRatio",naming2,sep=""),individualRatio)
  
 
  # incidence individual
   individualIncDiff <- riskdifference(Inc_11, #incidence (subtract initial HIV pop)
                 Inc_01, #incidence
                 N_11[HIV_11$t != 0],N_01[HIV_11$t != 0])
  tbl2[1,2] <- individualIncDiff$estimate # risk difference for incidence
  #confidence intervals
  assign(paste("lowerIncCI"), quantile(individualIncDiff$estimate, probs = .025))
  assign(paste("upperIncCI"),quantile(individualIncDiff$estimate, probs = .975))
  tbl2[1,3] <- paste(lowerIncCI,upperIncCI)
  indIncRatio <- riskratio(Inc_11, #incidence (subtract initial HIV pop)
                           Inc_01, #incidence
                           N_11,N_01)
  tbl2[1,4] <- ave(individualIncRatio$estimate)[1]
  tbl2[1,5] <- paste(quantile(individualIncRatio$estimate, probs = .025),
                     quantile(individualIncRatio$estimate, probs = .975))
  assign(paste("individualDiff",naming4,sep=""),indIncDiff)
  assign(paste("individualRatio",naming4,sep=""),indIncRatio)
  
# disseminated prevalence risk
  dissDiff <- riskdifference(HIV_01,HIV_00,
                             N_01,N_00)
  tbl1[2,2] <- dissDiff$estimate
  tbl1[2,3] <- paste(quantile(dissDiff$estimate, probs = .025),
                     quantile(dissDiff$estimate, probs = .975))
  
  dissRatio <- riskratio(HIV_01,HIV_00,
                         N_01,N_00)
  tbl1[2,4] <- dissRatio$estimate
  tbl1[2,5] <- paste(quantile(dissRatio$estimate, probs = .025),
                     quantile(dissRatio$estimate, probs = .975))
  
  assign(paste("disseminatedDiff",naming2, sep=""),dissDiff)
  assign(paste("disseminatedRatio",naming2, sep=""),dissRatio)
  
  #disseminated incident risk
  dissIncDiff <- riskdifference(Inc_01,
                                Inc_00,
                                N_01,N_00)
  tbl2[2,2] <- ave(dissIncDiff$estimate)[1]
  tbl2[2,3] <- paste(quantile(dissIncRatio$estimate, probs = .025),
                     quantile(dissIncRatio$estimate, probs = .975))
  
  dissIncRatio <- riskratio(Inc_01,
                            Inc_00,
                            N_01,N_00)
  tbl2[2,4] <- ave(dissIncRatio$estimate)[1]
  tbl2[2,5] <-   tbl1[2,2] <- dissDiff$estimate
  tbl1[2,3] <- paste(quantile(dissIncRatio$estimate, probs = .025),
                     quantile(dissIncRatio$estimate, probs = .975))
  
  assign(paste("disseminatedDiff",naming4, sep=""),dissIncDiff)
  assign(paste("disseminatedRatio",naming4, sep=""),dissIncRatio)

  # overall prevalence risk
  
  overallDiff <- riskdifference((HIV_11+HIV_01),HIV_00,
                                (N_11+N_01),N_00)
  overallRatio <- riskratio((HIV_11+HIV_01),HIV_00,
                            (N_11+N_01),N_00)
  
  assign(print("overallDiff",naming2, sep=""),overallDiff)
  assign(print("overallRatio",naming2, sep=""),overallRatio)
  
  #overall incidence risk
  overallIncDiff <- riskdifference((Inc_11+Inc_01),Inc_00,
                                   (N_11+N_01),N_00)
  overallIncRatio <- riskratio((Inc_11+Inc_01),Inc_00,
                            (N_11+N_01),N_00)
  
  assign(print("overallDiff",naming4, sep=""),overallIncDiff)
  assign(print("overallRatio",naming4, sep=""),overallIncRatio)
  
  # composite prevalence risk
  compositeDiff <- riskdifference(HIV_11,HIV_00,
                                  N_11,N_00)
  tbl1[3,2] <- ave(compositeDiff$estimate)[1]
  tbl1[3,3] <- paste(quantile(compositeDiff$estimate, probs = .025),
                     quantile(compositeDiff$estimate, probs = .975))
  compositeRatio <- riskratio(HIV_11,HIV_00,
                              N_11,N_00)
  tbl1[3,4] <- ave(compositeRatio$estimate)[1]
  tbl1[3,5] <- paste(quantile(compositeRatio$estimate, probs = .025),
                     quantile(compositeRatio$estimate, probs = .975))
  
  assign(print("compositeDiff",naming2, sep=""),compositeDiff)
  assign(print("overallRatio",naming2, sep=""),compositeRatio)
  
  # composite incidence risk
  compositeIncDiff <- riskdifference(Inc_11,Inc_00,
                                  N_11,N_00)
  compositeRatio <- riskratio(Inc_11,Inc_00,
                              N_11,N_00)
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
