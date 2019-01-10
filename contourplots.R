library("fmsb")
library("knitr")
library("dplyr")
#import files



contourValuesDisseminated <- as.data.frame(matrix(nrow=9,ncol=9))
coverageLevels <- seq(.1,.9, by = .1)
rownames(contourValuesDisseminated) <- coverageLevels
colnames(contourValuesDisseminated) <- coverageLevels
contourValuesComposite <- contourValuesDisseminated
contourValuesOverall <- contourValuesDisseminated

#risk ratios and risk differences for disseminated
table2 <- as.data.frame(matrix(nrow = 4, ncol = 5))
colnames(table2) <- c('Effect', 'RD', '95% CI','RR','95% CI')
table2$Effect <- c('Direct','Disseminated','Composite','Overall')




tab_df(df,title = "Contours", col.header = matrixRowNames, show.rownames = T, alternate.rows = T,file="test.doc")

#need to write total person-years on/off PrEP and the number w/HIV for both

#start with disseminated, then do composite

#create x by y matrix
#riskdifference(a and b)
#put risk difference at matrix["a","b"]
#put negative of risk difference at matrix["b","a"]

differenceAndRatioContours <- function(risk1,risk2,mtx1,mtx2,mtx3){ #make this work for disseminated and composite

  disseminatedDiff <- mean(risk1$risk01_inc) - mean(risk2$risk01_inc)
  
  compositeDiff <- mean(risk1$risk11_inc) - mean(risk2$risk01_inc)

  overallDiff <- mean(risk1$risk1_inc) - mean(risk2$risk1_inc)
  
  name1 <- deparse(substitute(risk1))
  name2 <- deparse(substitute(risk2))
  
  xLocation <- as.numeric(substr(name1, nchar(name1)-1, nchar(name1)-1))
  yLocation <- as.numeric(substr(name2, nchar(name2)-1, nchar(name2)-1))
  
  mtx1[xLocation,yLocation] <- disseminatedDiff
  mtx1[yLocation,xLocation] <- -disseminatedDiff

  mtx2[xLocation,yLocation] <- compositeDiff
  mtx2[yLocation,xLocation] <- -compositeDiff
  mtx3[xLocation,yLocation] <- overallDiff
  mtx3[yLocation,xLocation] <- -overallDiff
  #will this work??
  m1name <- deparse(substitute(mtx1))
  assign("contourValuesDisseminated", value = mtx1, envir = .GlobalEnv)
  assign("contourValuesComposite", value = mtx2, envir = .GlobalEnv)
  assign("contourValuesOverall", value = mtx3, envir = .GlobalEnv)
  contourValuesDisseminated <<- mtx1
}

differenceAndRatioContours(riskDF_10, riskDF_10, contourValuesDisseminated, contourValuesComposite, contourValuesOverall)
differenceAndRatioContours(riskDF_10, riskDF_20, contourValuesDisseminated, contourValuesComposite, contourValuesOverall)

differenceAndRatioContours(get(paste0("riskDF_",*,"0")), riskDF_20, contourValuesDisseminated, contourValuesComposite, contourValuesOverall)

riskNames <- list()
for (i in 1:9){
  riskNames[i] <- (paste0("riskDF_",as.character(i),"0"))
}
riskNames <- as.character(riskNames)
for (i in riskNames){
  for (j in riskNames){

    assign("r11", get(riskNames[1]))
    assign("r21", get(riskNames[2]))
    differenceAndRatioContours(r11,r21,contourValuesDisseminated, contourValuesComposite,
                               contourValuesOverall)
    contourValuesDisseminated <<- contourValuesDisseminated
    #assign("contourValuesDisseminated", value = contourValuesDisseminated, envir = .GlobalEnv)
    #assign("contourValuesComposite", value = contourValuesComposite, envir = .GlobalEnv)
    #assign("contourValuesOverall", value = contourValuesOverall, envir = .GlobalEnv)
  }
}
plot <- contour(x = rownames(contourValuesDisseminated),colnames(contourValuesDisseminated), contourValuesDisseminated)
#ggsave()