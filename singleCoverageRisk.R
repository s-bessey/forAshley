# This code will generate tables 2 and 4, plus figures 3 and 4. All comparisons 
# are to the control case (here no PrEP)

table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% CI','RR','95% CI')
table2$Effect <- c('Individual','Disseminated','Composite','Overall')


riskDiffRatio <- function(treatment){
  
  naming <- substr(treatment, nchar(treatment)-2,nchar(treatment))
  
  individualDiff <- riskdifference(treatment$HIV_PrEP,treatment$HIV_noPrEP,
                             treatment$PrEP,treatment$noPrEP)
  individualRatio <- riskratio(treatment$HIV_PrEP,treatment$HIV_noPrEP,
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
table2[1,2] <- individualDiff$estimate #rename all of these with proper % data
table2[1,3] <- individualDiff$conf.int
table2[1,4] <- individualRatio$estimate
table2[1,5] <- individualRatio$conf.int

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

