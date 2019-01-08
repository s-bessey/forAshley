library("fmsb")
library("knitr")
library("dplyr")
#import files



contourValuesDisseminated <- as.data.frame(matrix(nrow=9,ncol=9))
rownames(contourValuesDisseminated) <- seq(.1,.9,by=.1)
colnames(contourValuesDisseminated) <- seq(.1,.9,by=.1)
contourValuesCumulative <- contourValuesDisseminated
contourValuesOverall <- contourValuesDisseminated

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

differenceAndRatioContours <- function(risk1,risk2,mtx1,mtx2,mtx3){ #make this work for disseminated and composite

  disseminatedDiff <- mean(risk1$risk01_inc) - mean(risk2$risk00_inc)

  compositeDiff <- mean(risk1$risk11_inc) - mean(risk2$risk00_inc)
  
  overallDiff <- mean(risk1$risk1_inc) - mean(risk2$risk00_inc)
  name1 <- deparse(substitute(risk1))
  name2 <- deparse(substitute(risk2))
  xLocation <- substr(name1, nchar(name1)-1, nchar(name1)-1)
  #xLocation <- paste(".",xLocation) #may not need this if coverage is decimal?
  yLocation <- substr(name2, nchar(name2)-1, nchar(name2)-1)
  #yLocation <- paste(".",yLocation) #see above
  mtx1[xLocation,yLocation] <- disseminatedDiff
  #mtx1[yLocation,xLocation] <- -disseminatedDiff
  #mtx2[xLocation,yLocation] <- compositeDiff
  #mtx2[ylocation,xLocation] <- -compositeDiff
  #will this work??
  m1name <- deparse(substitute(mtx1))
  assign("contourValuesDisseminated", value = mtx1, envir = .GlobalEnv)
}

differenceAndRatioContours(riskDF_10, riskDF_20, contourValuesDisseminated, contourValuesCumulative, contourValuesOverall)

for (i in 1:9){
  for (j in 1:9){
    differenceAndRatioContours(paste("riskDF_",i,"0", sep = ""),paste("riskDF_",j,"0", sep = ""),
                               contourValuesDisseminated, contourValuesCumulative,
                               contourValuesOverall)
  }
}
plot = contour(x = matrixRowNames,y = matrixRowNames, contourValues)
#ggsave()