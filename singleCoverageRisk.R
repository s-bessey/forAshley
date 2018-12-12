# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)


#have to have cum and prev ave first
table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% CI','RR','95% CI')
table2$Effect <- c('Individual','Disseminated','Composite','Overall')
table4 <- table2

riskDiffRatioT2 <- function(treatment,tbl1,tbl2){ #give the file of endtime data
  # and prevalence or incidence
  
  naming2 <- print(substr(treatment, nchar(treatment)-2,nchar(treatment)),
                  "Prev", sep = "")
  naming4 <- print(substr(treatment, nchar(treatment)-2,nchar(treatment)),
                   "Inc", sep = "")
  
  N_11 <- treatment$PrEP # number of people on PrEP
  HIV_11 <- treatment$N_HIV_PrEP # prevalence of HIV+ PrEP-users
  Inc_11 <- treatment$Inc_HIV_PrEP # incidence of HIV in PrEP users
    
  N_01 <- treatment$N_noPrEP # number of people not on PrEP in treatment gp
  HIV_01 <- treatment$HIV_noPreP # HIV+ people not on PrEP in treatment gp
  Inc_01 <- treatment$Inc_HIV_noPrEP # HIV incidence in no-PrEP treatment gp
  
  N_00 <- treatment$N_control # population of control gp
  HIV_00 <- treatment$HIV_control # HIV prevalence in control
  Inc_00 <- treatment$Inc_control # HIV incidence in control
  
  
  # create differences and ratios for both prevalence and incidence,
  # 
  individualDiff <- riskdifference(HIV_11,HIV_01,
                             N_11,N_01)
  tbl1[1,2] <- individualDiff$estimate
  tbl1[1,3] <- paste(individualDiff$conf.int[1],individualDiff$conf.int[2])
  individualRatio <- riskratio(HIV_11,HIV_01,
                               N_11,N_01)
  tbl1[1,4] <- individualRatio$estimate
  tbl1[1,5] <- paste(individualRatio$conf.int[1],individualRatio$conf.int[2])
  assign(paste("individualDiff",naming2,sep=""),individualDiff)
  assign(paste("individualRatio",naming2,sep=""),individualRatio)
  
  indIncDiff <- riskdifference(Inc_11,Inc_01,
                                N_11,N_01)
  tbl2[1,2] <- individualIncDiff$estimate
  tbl2[1,3] <- paste(individualIncDiff$conf.int[1],individualIncDiff$conf.int[2])
  indIncRatio <- riskratio(Inc_11,Inc_01,
                                N_11,N_01)
  tbl2[1,4] <- individualIncRatio$estimate
  tbl2[1,5] <- paste(individualIncRatio$conf.int[1],individualIncRatio$conf.int[2])
  assign(paste("individualDiff",naming4,sep=""),indIncDiff)
  assign(paste("individualRatio",naming4,sep=""),indIncRatio)
  

  dissDiff <- riskdifference(HIV_01,HIV_00,
                             N_01,N_00)
  tbl1[2,2] <- dissDiff$estimate
  tbl1[2,3] <- paste(dissDiff$conf.int[1],dissDiff$conf.int[2])
  dissRatio <- riskratio(HIV_01,HIV_00,
                         N_01,N_00)
  tbl1[2,4] <- dissRatio$estimate
  tbl1[2,5] <- paste(dissRatio$conf.int[1],dissRatio$conf.int[2])
  
  assign(paste("disseminatedDiff",naming2, sep=""),dissDiff)
  assign(paste("disseminatedRatio",naming2, sep=""),dissRatio)
  
  dissIncDiff <- riskdifference(Inc_01,Inc_00,
                             N_01,N_00)
  tbl2[2,2] <- dissIncDiff$estimate
  tbl2[2,3] <- paste(dissIncDiff$conf.int[1],dissIncDiff$conf.int[2])
  dissIncRatio <- riskratio(Inc_01,Inc_00,
                         N_01,N_00)
  tbl2[2,4] <- dissIncRatio$estimate
  tbl2[2,5] <- paste(dissIncRatio$conf.int[1],dissIncRatio$conf.int[2])
  
  assign(paste("disseminatedDiff",naming4, sep=""),dissIncDiff)
  assign(paste("disseminatedRatio",naming4, sep=""),dissIncRatio)

  # everything under this needs correction
  
  overallDiff <- riskdifference((HIV_11+HIV_01),HIV_00,
                                (N_11+N_01),N_00)
  overallRatio <- riskratio((HIV_11+HIV_01),HIV_00,
                            (N_11+N_01),N_00)
  
  assign(print("overallDiff",naming2, sep=""),overallDiff)
  assign(print("overallRatio",naming2, sep=""),overallRatio)
  
  overallIncDiff <- riskdifference((Inc_11+Inc_01),Inc_00,
                                   (N_11+N_01),N_00)
  overallIncRatio <- riskratio((Inc_11+Inc_01),Inc_00,
                            (N_11+N_01),N_00)
  
  assign(print("overallDiff",naming4, sep=""),overallIncDiff)
  assign(print("overallRatio",naming4, sep=""),overallIncRatio)
  
  compositeDiff <- riskdifference(HIV_11,HIV_00,
                                  N_11,N_00)
  compositeRatio <- riskratio(HIV_11,HIV_00,
                              N_11,N_00)
  
  assign(print("compositeDiff",naming2, sep=""),compositeDiff)
  assign(print("overallRatio",naming2, sep=""),compositeRatio)
  
  compositeIncDiff <- riskdifference(Inc_11,Inc_00,
                                  N_11,N_00)
  compositeRatio <- riskratio(Inc_11,Inc_00,
                              N_11,N_00)
  assign(print("compositeDiff",naming4, sep=""),compositeIncDiff)
  assign(print("overallRatio",naming4, sep=""),compositeIncRatio)
}

riskDiffRatio(aveTreatment10)
#make data frame for prevalence
table2[1,2] <- individualDiff$estimate #rename all of these with proper % data
table2[1,3] <- paste(individualDiff$conf.int[1],individualDiff$conf.int[2])
table2[1,4] <- individualRatio$estimate
table2[1,5] <- paste(individualRatio$conf.int[1],individualRatio$conf.int[2])

table2[2,2] <- disseminatedDiff$estimate
table2[2,3] <- disseminatedDiff$conf.int
table2[2,4] <- disseminatedRatio$estimate
table2[2,5] <- disseminatedRatio$conf.int

table2[3,2] <- compositeDiff$estimate
table2[3,3] <- compositeDiff$conf.int
table2[3,4] <- compositeRatio$estimate
table2[3,5] <- compositeRatio$conf.int

table2[4,2] <- overallDiff$estimate
table2[4,3] <- overallDiff$conf.int
table2[4,4] <- overallRatio$estimate
table2[4,5] <- overallRatio$conf.int

tab_df(table2, col.header = colnames(table2), show.rownames = F, alternate.rows = T,file="Table2.doc")

ggplot(individualEffects, aes(x = level, y= estimate)) +
  geom_ribbon(aes(ymin = CILower,ymax = CIUpper))

