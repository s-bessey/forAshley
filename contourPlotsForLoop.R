differenceAndRatioContoursTest <- function(mtx1,mtx2,mtx3){ #make this work for disseminated and composite
  riskNames <- list()
  for (i in 1:9){
  riskNames[i] <- (paste0("riskDF_",as.character(i),"0"))
  }
  for (i in 1:9){

    for (j in 1:9){
      assign("r11", get(as.character(riskNames[i])))
      #r21 <- get(riskNames[j])
      #differenceAndRatioContours(r11,r21,contourValuesDisseminated, contourValuesComposite,
      #                           contourValuesOverall)
      disseminatedDiff <- mean(r11$risk01_inc) - mean(r21$risk01_inc)
      name1 <- riskNames[i]
      name2 <- riskNames[j]
      xLocation <- as.numeric(substr(name1, nchar(name1)-1, nchar(name1)-1))
      yLocation <- as.numeric(substr(name2, nchar(name2)-1, nchar(name2)-1))
      
      mtx1[xLocation,yLocation] <- disseminatedDiff
      mtx1[yLocation,xLocation] <- -disseminatedDiff
      assign("contourValuesDisseminated", value = mtx1, envir = .GlobalEnv)
    }
      
  }
  
  
return(name1)
  #disseminatedDiff <- mean(risk1$risk01_inc) - mean(risk2$risk01_inc)
  
  xLocation <- as.numeric(substr(name1, nchar(name1)-1, nchar(name1)-1))
  yLocation <- as.numeric(substr(name2, nchar(name2)-1, nchar(name2)-1))
  
  mtx1[xLocation,yLocation] <- disseminatedDiff
  mtx1[yLocation,xLocation] <- -disseminatedDiff
  
  #mtx2[xLocation,yLocation] <- compositeDiff
  #mtx2[yLocation,xLocation] <- -compositeDiff
  #mtx3[xLocation,yLocation] <- overallDiff
  #mtx3[yLocation,xLocation] <- -overallDiff
  #will this work??
  #m1name <- deparse(substitute(mtx1))
  assign("contourValuesDisseminated", value = mtx1, envir = .GlobalEnv)
  #assign("contourValuesComposite", value = mtx2, envir = .GlobalEnv)
  #assign("contourValuesOverall", value = mtx3, envir = .GlobalEnv)
  #contourValuesDisseminated <<- mtx1
}
