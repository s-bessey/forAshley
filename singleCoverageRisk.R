# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)


#have to have cum and prev ave first
table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% CI','RR','95% CI')
table2$Effect <- c('Individual','Disseminated','Composite','Overall')
table4 <- table2

riskDiffRatioT2 <- function(treatment){ #give the file of endtime data
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
  
  
  
  individualDiff <- riskdifference(treatment$HIV_PrEP_prev,treatment$HIV_noPrEP_prev,
                             treatment$PrEP_prev,treatment$noPrEP_prev)
  individualRatio <- riskratio(treatment$HIV_PrEP_prev,treatment$HIV_noPrEP_prev,
                               treatment$PrEP,treatment$noPrEP)
  assign(print("individulaDiff",naming,sep=""),individualDiff)
  assign(print("individulaRatui",naming,sep=""),individualRatio)
  
  dissDiff <- riskdifference(treatment$HIV_noPrEP,treatment$HIV_control,
                             treatment$noPrEP,treatment$control)
  dissRatio <- riskratio(treatment$HIV_noPrEP,treatment$HIV_control,
                         treatment$noPrEP,treatment$control)
  
  assign(print("disseminatedDiff",naming, sep=""),dissDiff)
  assign(print("disseminatedRatio",naming, sep=""),dissRatio)
  
  overallDiff <- riskdifference(treatment$HIV,treatment$HIV_control,
                                (treatment$PrEP+treatment$noPrEP),treatment$control)
  overallRatio <- riskratio(treatment$HIV,treatment$HIV_control,
                            (treatment$PrEP+treatment$noPrEP),treatment$control)
  assign(print("overallDiff",naming, sep=""),overallDiff)
  assign(print("overallRatio",naming, sep=""),overallRatio)
  
  compositeDiff <- riskdifference(treatment$HIV_PrEP,treatment$HIV_control,
                                  treatment$PrEP,treatment$control)
  compositeRatio <- riskratio(treatment$HIV_PrEP,treatment$HIV_control,
                              treatment$PrEP,treatment$control)
  assign(print("compositeDiff",naming, sep=""),compositeDiff)
  assign(print("overallRatio",naming, sep=""),compositeRatio)
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

