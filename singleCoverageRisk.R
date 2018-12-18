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

riskDiffRatio <- function(treatment,tbl1,tbl2,coverage){ #give the file of endtime data
  # and prevalence or incidence
  #may move the risks down when using estimators
  naming2 <- paste(coverage,"Prev", sep = "") # name for prevalence
  naming4 <- paste(coverage,"Inc", sep = "") # name for incidence
 
  N_11 <- treatment$Nprep[treatment$t == 0 & treatment$NprepElig !=0] # number of people on PrEP
  HIV_11 <- integer(length = length(N_11)) # no PrEP users have HIV
  Prev_risk_11 <- HIV_11 / N_11
  Prev_11_quants <- quantile(Prev_risk_11,probs = .025,.95)
  Inc_11 <- integer(length = length(N_11))
  Inc_risk_11 <- Inc_11 / N_11
  Inc_11_quants <- quantile(Inc_risk_11,probs = .025,.95)

  time0 <- subset(treatment,t == 0)
  time60 <- subset(treatment, t == max(treatment$t))
  start_01 <- subset(time0, NprepElig != 0) # number of people not on PrEP in treatment gp
  N_01 <- start_01$NprepElig
  end_01 <- subset(time60, time0$NprepElig != 0)
  HIV_01 <- end_01$Nhiv

  Prev_risk_01 <- HIV_01 / N_01
  Prev_01_quants <- quantile(Prev_risk_01,probs = .025,.95)
  
  Inc_01 <- end_01$Nhiv # HIV incidence in people in treatment gp with no PrEP

  Inc_risk_01 <- Inc_01 / N_01
  Inc_01_quants <- quantile(Inc_risk_01,probs = .025,.95)


  N_00 <- treatment$N[treatment$NprepElig == 0 & treatment$Nprep == 0 & treatment$t == 0] # population of control gp
  HIV_00 <- treatment$Nhiv[treatment$NprepElig == 0 & treatment$Nprep == 0 & treatment$t == 0] # HIV population of control gp

  Prev_risk_00 <- HIV_00 / N_00
  Prev_00_quants <- quantile(Prev_risk_00,probs = .025,.95)
  Inc_00 <- treatment$Nnewinf[treatment$t == max(treatment$t) & treatment$NprepElig == 0 &
                                    treatment$Nprep == 0] # HIV incidence in control
  Inc_risk_00 <- Inc_00 / N_00
  Inc_00_quants <- quantile(Inc_risk_00,probs = .025,.95)
  

  # create differences and ratios for both prevalence and incidence
  # individual prevalence
  individualDiff <- ave(Prev_risk_11)[1] - ave(Prev_risk_01)[1]
  tbl1[1,2] <- individualDiff #add to table
  #lowerPrevCI
  lowerPrevCI<-Prev_11_quants[1]-Prev_01_quants[2]
  upperPrevCI <- Prev_11_quants[2]-Prev_01_quants[1]
  tbl1[1,3] <- paste(lowerPrevCI,upperPrevCI) #find simulation intervals (2.5 and 97.5 percentiles)
  individualRatio <- ave(Prev_risk_11)[1]/ave(Prev_risk_01)[1]

  tbl1[1,4] <- individualRatio
  lowerIncCI <- Prev_11_quants[1]/Prev_01_quants[2]
  upperIncCI <- Prev_11_quants[2]/Prev_01_quants[1]
  tbl1[1,5] <- paste(lowerIncCI,upperIncCI)
  assign(paste("individualDiff",naming2,sep=""),individualDiff)
  assign(paste("individualRatio",naming2,sep=""),individualRatio)


  # incidence individual
  individualIncDiff <- ave(Inc_risk_11)[1]-ave(Inc_risk_01)
  tbl2[1,2] <- individualIncDiff # risk difference for incidence
  #confidence intervals
  lowerIncCI <- Inc_11_quants[1]-Inc_01_quants[2]
  upperIncCI <- Inc_11_quants[2]-Inc_01_quants[1]
  tbl2[1,3] <- paste(lowerIncCI,upperIncCI)
  indIncRatio <- ave(Inc_risk_11)[1]/ave(Inc_risk_01)[1]
  lowerIncCI <- Inc_11_quants[1]/Inc_01_quants[2]
  upperIncCI <- Inc_11_quants[2]/Inc_01_quants[1]
  tbl2[1,4] <- individualIncRatio
  tbl2[1,5] <- paste(lowerIncCI,upperIncCI)
  assign(paste("individualDiff",naming4,sep=""),indIncDiff)
  assign(paste("individualRatio",naming4,sep=""),indIncRatio)
  
  
  #-----
  #fix difference and ratios beyond this point
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
