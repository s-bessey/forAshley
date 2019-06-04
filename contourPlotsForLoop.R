contourValuesDisseminated <- as.data.frame(matrix(nrow=9,ncol=9))
coverageLevels <- seq(.1,.9, by = .1)
rownames(contourValuesDisseminated) <- coverageLevels
colnames(contourValuesDisseminated) <- coverageLevels
contourValuesComposite <- contourValuesDisseminated
contourValuesOverall <- contourValuesDisseminated

differenceAndRatioContoursTest <- function(mtx1,mtx2,mtx3){
  riskNames <- list()
  for (i in 1:9){
  riskNames[i] <- (paste0("risk",as.character(i),"0"))
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
      
      mtx1[xLocation,yLocation] <- -disseminatedDiff
      mtx1[yLocation,xLocation] <- disseminatedDiff
      mtx2[xLocation,yLocation] <- -compositeDiff
      mtx2[yLocation,xLocation] <- compositeDiff
      mtx3[xLocation,yLocation] <- -overallDiff
      mtx3[yLocation,xLocation] <- overallDiff
      assign("contourValuesDisseminated", value = mtx1, envir = .GlobalEnv)
      assign("contourValuesComposite", value = mtx2, envir = .GlobalEnv)
      assign("contourValuesOverall", value = mtx3, envir = .GlobalEnv)
    }
      
  }
}
differenceAndRatioContoursTest(contourValuesDisseminated,contourValuesComposite,contourValuesOverall)

cols = rev(colorRampPalette(c("slateblue4",'slategray1'))(20))
#cols = rev(colorRampPalette(c('#e66101', '#fdb863', '#b2abd2', '#5e3c99'))(20))

png("DisseminatedContourFilled0603.png",res=600,height=8.5,width=11,units="in")
filled.contour(x = coverageLevels, y = coverageLevels, z = as.matrix(contourValuesDisseminated),
               plot.axes = {contour(coverageLevels,coverageLevels,as.matrix(contourValuesDisseminated),nlevels=14,drawlabels=TRUE, axes=TRUE,add=TRUE); axis(side = 1,at=seq(.2,.8,by=.1));axis(side=2,at=seq(.2,.8,by=.1))},
               xlab =expression(alpha), ylab = expression(paste(alpha,"'",sep = "")),axes=T)
#contour(x = coverageLevels, y = coverageLevels, z = as.matrix(contourValuesDisseminated),
#        xlab =expression(alpha), ylab = expression(paste(alpha,"'",sep = "")))
dev.off()


png("OverallContour0603.png",res=600,height=8.5,width=11,units="in")
filled.contour(x = coverageLevels, y = coverageLevels,
        z = as.matrix(contourValuesOverall), nlevels=14, #col = cols,
        plot.axes = {contour(coverageLevels,coverageLevels,as.matrix(contourValuesOverall),nlevels=14,drawlables=TRUE, axes=TRUE,add=TRUE); axis(side = 1,at=seq(.2,.8,by=.1));axis(side=2,at=seq(.2,.8,by=.1))},
        xlab =expression(alpha),
        ylab = expression(paste(alpha,"'",sep = "")))
dev.off()

png("CompositeContourFilled0603.png",res=600,height=8.5,width=11,units="in")
filled.contour(x = coverageLevels, y = coverageLevels,
        z = as.matrix(contourValuesComposite),nlevels=14, #col = cols,
        plot.axes = {contour(coverageLevels,coverageLevels,as.matrix(contourValuesComposite),nlevels=14,drawlables=TRUE, axes=TRUE,add=TRUE); axis(side = 1,at=seq(.2,.8,by=.1));axis(side=2,at=seq(.2,.8,by=.1))},
        xlab =expression(alpha),
        ylab = expression(paste(alpha,"'",sep = "")))

dev.off()

