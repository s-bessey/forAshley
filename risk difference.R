library("fmsb")
library("knitr")
library("dplyr")
#import files

dfTreatment10 <- read.table() #insert filename

aveTreatment10 <- aggregate(.~t, dfTreatment10,function(x) mean = mean(x))

matrixRowNames <- seq(.1,.9,by=.1)
contourValuesDisseminated <- matrix(nrow=9,ncol=9)
contourValuesCumulative <- contourValuesDisseminated
rownames(contourValues) <- matrixRowNames
colnames(contourValues) <- matrixRowNames

#risk ratios and risk differences for disseminated
table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% CI','RR','95% CI')
table2$Effect <- c('Direct','Disseminated','Composite','Overall')




riskDiffRatio(aveTreatment10)

tab_df(df,title = "Contours", col.header = matrixRowNames, show.rownames = T, alternate.rows = T,file="test.doc")

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