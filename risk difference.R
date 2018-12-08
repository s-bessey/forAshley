library("fmsb")
#import files

dfTreatment10 <- read.table() #insert filename

aveTreatment10 <- aggregate(.~t, dfTreatment10,function(x) mean = mean(x))

matrixRowNames <- seq(.1,.9,by=.1)
contourValues <- matrix(nrow=9,ncol=9)
rownames(contourValues) <- matrixRowNames
colnames(contourValues) <- matrixRowNames

RD10 <- riskdifference(aveTreatment10$HIV_PrEP,aveTreatment10$HIV_PrEP,
                       aveTreatment10$PrEP, aveTreatment10$noPrEP)

RD10 <- riskratio(aveTreatment10$HIV_PrEP,aveTreatment10$HIV_PrEP,
                       aveTreatment10$PrEP, aveTreatment10$noPrEP)
#need to write total person-years on/off PrEP and the number w/HIV for both

#start with disseminated, then do composite

#create x by y matrix
#riskdifference(a and b)
#put risk difference at matrix["a","b"]
#put negative of risk difference at matrix["b","a"]

incifferenceAndRatio <- function(treatment1,treatment2,mtx){ #make this work for disseminated and composite
  avg1 <- aggregate(.~t, treatment1,function(x) mean = mean(x))
  disseminatedDiff <- riskdifference(treatment1$HIV_noPrEP,treatment2$HIV_noPrEP,
                                     treatment1$noPrEP,treatment2$noPrEP)
  xLocation <- substr(treatment1, nchar(treatment1)-2, nchar(treatment1)-1)
  xLocation <- paste(".",xLocation)
  yLocation <- substr(treatment2, nchar(treatment2)-2, nchar(treatment2))
  yLocation <- paste(".",yLocation)
  mtx[xLocation,yLocation] <- disseminatedDiff
  #jesus christ will this work??
}