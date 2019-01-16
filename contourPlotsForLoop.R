differenceAndRatioContoursTest <- function(mtx1,mtx2,mtx3){
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
}
cols = rev(colorRampPalette(c("slateblue4",'slategray1'))(20))
cols = rev(colorRampPalette(c('#e66101', '#fdb863', '#b2abd2', '#5e3c99'))(20))

png("DisseminatedContour.png",res=600,height=8.5,width=11,units="in")
dissContour <- filled.contour(x = coverageLevels, y = coverageLevels, z = as.matrix(contourValuesDisseminated), col = cols,
                              xlab =expression(alpha), ylab = expression(paste(alpha,"'",sep = "")))
dev.off()


png("OverallContour.png",res=600,height=8.5,width=11,units="in")
overallContour <- filled.contour(x = coverageLevels, y = coverageLevels,
                                 z = as.matrix(contourValuesOverall), col = cols,
                                 xlab =expression(alpha),
                                 ylab = expression(paste(alpha,"'",sep = "")))
dev.off()

png("CompositeContour.png")
compositeContour <- filled.contour(x = coverageLevels, y = coverageLevels,
                                   z = as.matrix(contourValuesComposite), col = cols,
                                  xlab =expression(alpha),
                                  ylab = expression(paste(alpha,"'",sep = "")))
dev.off()

