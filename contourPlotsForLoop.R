differenceAndRatioContoursTest <- function(mtx1,mtx2,mtx3){ #make this work for disseminated and composite
  riskNames <- list()
  for (i in 1:9){
  riskNames[i] <- (paste0("riskDF_",as.character(i),"0"))
  }
  for (i in 1:9){

    for (j in 1:9){
      assign("r11", get(as.character(riskNames[i])))
      assign("r21", get(as.character(riskNames[j])))
      #r21 <- get(riskNames[j])
      #differenceAndRatioContours(r11,r21,contourValuesDisseminated, contourValuesComposite,
      #                           contourValuesOverall)
      disseminatedDiff <- mean(r11$risk01_inc) - mean(r21$risk01_inc)
      
      compositeDiff <- mean(r11$risk11_inc) - mean(r21$risk01_inc)
      
      overallDiff <- mean(r11$risk1_inc) - mean(r21$risk1_inc)
      
      name1 <- riskNames[i]
      name2 <- riskNames[j]
      xLocation <- as.numeric(substr(name1, nchar(name1)-1, nchar(name1)-1))
      yLocation <- as.numeric(substr(name2, nchar(name2)-1, nchar(name2)-1))
      
      mtx1[xLocation,yLocation] <- disseminatedDiff
      mtx1[yLocation,xLocation] <- -disseminatedDiff
      mtx2[xLocation,yLocation] <- compositeDiff
      mtx2[yLocation,xLocation] <- -compositeDiff
      mtx3[xLocation,yLocation] <- overallDiff
      mtx3[yLocation,xLocation] <- -overallDiff
      assign("contourValuesDisseminated", value = mtx1, envir = .GlobalEnv)
      assign("contourValuesComposite", value = mtx2, envir = .GlobalEnv)
      assign("contourValuesOverall", value = mtx3, envir = .GlobalEnv)
    }
      
  }


dissContour <- filled.contour(x = coverageLevels, y = coverageLevels, z = as.matrix(contourValuesDisseminated))
ggsave("DisseminatedContour.pdf", dissContour)
overallContour <- filled.contour(x = coverageLevels, y = coverageLevels, z = as.matrix(contourValuesOverall))
ggsave("OverallContour.pdf",overallContour)
compositeContour <- filled.contour(x = coverageLevels, y = coverageLevels, z = as.matrix(contourValuesComposite))
ggsave("CompositeContour.pdf", compositeContour)
