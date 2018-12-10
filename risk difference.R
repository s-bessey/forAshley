library("fmsb")
#import files

dfTreatment10 <- read.table() #insert filename

aveTreatment10 <- aggregate(.~t, dfTreatment10,function(x) mean = mean(x))

matrixRowNames <- seq(.1,.9,by=.1)
contourValuesDisseminated <- matrix(nrow=9,ncol=9)
contourValuesCumulative <- contourValuesDisseminated
rownames(contourValues) <- matrixRowNames
colnames(contourValues) <- matrixRowNames

#risk ratios and risk differences for disseminated

riskDiffRatio <- function(treatment){
  
  naming <- substr(treatment, nchar(treatment)-2,nchar(treatment))
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


#need to write total person-years on/off PrEP and the number w/HIV for both

#start with disseminated, then do composite

#create x by y matrix
#riskdifference(a and b)
#put risk difference at matrix["a","b"]
#put negative of risk difference at matrix["b","a"]

differenceAndRatioContours <- function(treatment1,treatment2,mtx1,mtx2){ #make this work for disseminated and composite
  avg1 <- aggregate(.~t, treatment1,function(x) mean = mean(x))
  disseminatedDiff <- riskdifference(treatment1$HIV_noPrEP,treatment2$HIV_noPrEP,
                                     treatment1$noPrEP,treatment2$noPrEP)
  compositeDiff <- riskdifference(treatment$HIV_PrEP,treatment$HIV_control,
                                  treatment$PrEP,treatment$control)
  xLocation <- substr(treatment1, nchar(treatment1)-2, nchar(treatment1)-1)
  xLocation <- paste(".",xLocation) #may not need this if coverage is decimal?
  yLocation <- substr(treatment2, nchar(treatment2)-2, nchar(treatment2))
  yLocation <- paste(".",yLocation) #see above
  mtx1[xLocation,yLocation] <- disseminatedDiff
  mtx1[yLocation,xLocation] <- -disseminatedDiff
  mtx2[xLocation,yLocation] <- compositeDiff
  mtx2[ylocation,xLocation] <- -compositeDiff
  #will this work??
}

plot = contour(x = matrixRowNames,y = matrixRowNames, contourValues)
#ggsave()