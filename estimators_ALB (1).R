N<-11245
# total number of agents in run isn't in the file, so set it here for ease
# of changing as necessary
library(dplyr)
treatment10<- read.table("reduced_24_0.1.txt",header = T)
treatment20<- read.table("reduced_24_0.2.txt",header = T)
treatment30<- read.table("reduced_24_0.3.txt",header = T)
treatment40<- read.table("reduced_24_0.4.txt",header = T)
treatment50<- read.table("reduced_24_0.5.txt",header = T)
treatment60<- read.table("reduced_24_0.6.txt",header = T)
treatment70<- read.table("reduced_24_0.7.txt",header = T)
treatment80<- read.table("reduced_24_0.8.txt",header = T)
treatment90<- read.table("reduced_24_0.9.txt",header = T)

for (i in 1:ncol(treatment)){
  for (j in 1:nrow(treatment)){
    if(is.character(treatment[i,j])){
      treatment10[i,j]<-as.numeric(treatment[i,j])
    }
  }
}
#setwd("/Volumes/GoogleDrive/My Drive/Networks/Brandon/Model Results/Ashley_runs_Dec2018")
temp=list.files(pattern = "*10_*")
# This should be a function, but it was easier to do this way
# these read the names of everything for a given coverage level, then bind those
# into one data frame
datalist <- list()
temp=list.files(pattern = "*10_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  temp2 <- subset(temp2, select = -c(pseed, rseed))
  temp2 <- temp2[which(temp2$t == max(temp2$t)|temp2$t == 0),]
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment10 <- do.call(rbind, datalist)


temp=list.files(pattern = "*20_*")
for (i in 1:10){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE, sep = "", stringsAsFactors = F))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment20 <-do.call(rbind,datalist)
datalist <- list()

temp=list.files(pattern = "*30_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist <- NULL
  datalist[[i]] <- temp2
}
treatment30 <-do.call(rbind,datalist)
datalist <- NULL

temp=list.files(pattern = "*40_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment40 <-do.call(rbind,datalist)
datalist <-NULL

temp=list.files(pattern = "Report_50_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment50 <-do.call(rbind,datalist)

temp=list.files(pattern = "*60_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment60 <-do.call(rbind,datalist)


temp=list.files(pattern = "*70_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment70 <-do.call(rbind,datalist)


temp=list.files(pattern = "*80_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment80 <-do.call(rbind,datalist)


temp=list.files(pattern = "*90_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment90 <-do.call(rbind,datalist)


#K = 2830,2822,1468,1438,1490,1495,1517,1486, 1489
####FUNCTION#####
#this calculates the parts used in the estimators
estimators <- function(treatment,coverage){
  # adds a column of zeroes (non-treatment) or ones (treatment)
  # should have used an if statement, can be fixed later

  #treatment$treatmentgp <- (treatment$Nprep + treatment$NprepElig)/(treatment$Nprep + treatment$NprepElig)
  #treatment[is.na(treatment)] <-0
  
  treatment$treatmentgp<-ifelse(treatment$NprepElig>0,1,0)
  
  #subset based on time so that things aren't calculated at t= 0 AND t = 60
  time0 <- subset(treatment,t == 0)
  time60 <- subset(treatment, t == max(treatment$t))
  # main data frame starts with the seed and component ID
  infectionValues <- as.data.frame(cbind(time60$nseed,time60$compID))
  colnames(infectionValues) <- c("nseed","compID")
  
  # start to append needed values for calculations to the data frame
  # find total people in each component (t = 0)
  infectionValues$totalN <- time0$totalN
  #denominator for prevalence
  # lose prep elig = no. prep elig infected
  infectionValues$N_11 <- time0$Nprep*time0$treatmentgp # number of people on PrEP
  # number of non-treated people in treatment group at time 0
  # includes PrEP-eligible and those who start infected
  infectionValues$N_01 <- (time0$totalN-time0$Nprep)*time0$treatmentgp
  
  #denominator for incidence (only among those not HIV+ at time =0)
  infectionValues$N_11i <- time0$Nprep*time0$treatmentgp
  infectionValues$N_01i <- (time0$NprepElig) #why am i subtracting number on prep???
  infectionValues$totalNi <- time0$totalN-time0$Nhiv
  # people lose PrEP eligibility by becoming infected. Other infections are
  # those who started the model on PrEP. So # who stop being eligible + number
  # who started HIV+ = prevalence among untreated treatment group
  # infectionValues$HIV_01 <- (time0$NprepElig - time60$NprepElig) + (time0$Nhiv*time0$treatmentgp)
  # # see above
  # infectionValues$Inc_01 <- time0$NprepElig - time60$NprepElig # no. of people who
  # # see above re: treatment group infections
  # infectionValues$HIV_11 <- 0
  # infectionValues$Inc_11 <- (time60$Nnewinf * time0$treatmentgp) - infectionValues$Inc_01

  
  infectionValues$HIV_11 <- time0$treatmentgp*(time60$NprepeverHIV)
  infectionValues$HIV_01 <- time0$treatmentgp*(time60$Nhiv)-infectionValues$HIV_11
  infectionValues$Inc_11 <- time0$treatmentgp*(time60$NprepeverHIV) #edit
  infectionValues$Inc_01 <- time0$treatmentgp*(time60$Nnewinf-infectionValues$Inc_11)
  
  infectionValues$N_trt_group <- infectionValues$N_11 + infectionValues$N_01
  infectionValues$N_trt_groupi <- infectionValues$N_11i + infectionValues$N_01i
  
  infectionValues$N_cntrl_group <- (time0$totalN)*(1-time0$treatmentgp)
  infectionValues$N_cntrl_groupi <- (time0$totalN-time0$Nhiv)*(1-time0$treatmentgp)
  
  # infectionValues$N_00 <- time0$totalN * (1- time0$treatmentgp) # population of control gp
  infectionValues$HIV_00 <- (time60$Nhiv)*(1- time0$treatmentgp) # HIV population of control gp
  infectionValues$Inc_00 <- (time60$Nhiv-time0$Nhiv)*(1- time0$treatmentgp) #Incidence in cntrl gp

  #individual weights
  # number of components (length of each nseed)
  # the average just makes it populate to everything in same nseed
  K <- with(time0,ave(compID,nseed,FUN = length))
  
  infectionValues$w_h <-  1/(infectionValues$totalN*K)
  infectionValues$w_i <- 1/N
 


  #find sums for each group
  #fist define each equation that's summed in the weight calcs
  #weights for incidence (among HIV-neg at time 0)
  # weightsi <- as.data.frame(cbind(infectionValues$nseed,infectionValues$N_trt_groupi,
  #                                infectionValues$N_11i*infectionValues$N_trt_groupi,
  #                             (1-infectionValues$N_11i)*infectionValues$N_trt_groupi,
  #                             (1-infectionValues$N_trt_groupi)*(infectionValues$N_trt_groupi))) 
  
  weightsi <- as.data.frame(cbind(infectionValues$nseed,
                                 infectionValues$N_trt_groupi,
                                 infectionValues$N_11i,
                                 infectionValues$N_01i,
                                 infectionValues$N_cntrl_groupi))
  colnames(weightsi) <- c("nseed", "sum1", "sum2","sum3","sum4")
  # take the sums needed
  sums1i <- aggregate(.~nseed,weightsi, function(x) sum = sum(x))
  # repopulate with a sums column the same length as the df
  sumsi <- merge(sums1i, infectionValues, by = "nseed")
  
  #infectionValues$sum1 <- tapply(infectionValues$N_trt_group,INDEX = infectionValues$nseed, sum)
  

  
  # two stage inverse probability weights
  # use sums from above
  sumsi$w11i <- (sumsi$sum1/N)^(-1) *(sumsi$sum2/sumsi$sum1)^(-1)
  #assign("test", w11i,envir = .GlobalEnv)

  sumsi$w01i <- (sumsi$sum1/N)^(-1) * (sumsi$sum3/sumsi$sum1)^(-1)
  sumsi$w0i <- (sumsi$sum4/N)^(-1)
  sumsi$w1i <- (sumsi$sum1/N)^(-1)
   #make sure the weights are in the correct order for data set
  sums2i<-as.data.frame(cbind(sumsi$nseed,sumsi$w11i,sumsi$w01i,sumsi$w0i,sumsi$w1i))
  names(sums2i) <- c("nseed", "w11","w01","w0","w1")
  
  #keep one record per nseed
  sums3i<- sums2i %>% 
    group_by(nseed) %>% 
    filter(row_number()==1) %>%
    as.data.frame()
  #sums2a <- do.call(rbind, lapply(split(sums2, sums2$nseed), head, 1))
  
  #weights for prevalence (among everyone at time 0)
  weights <- as.data.frame(cbind(infectionValues$nseed,
                                 infectionValues$N_trt_group,
                                 infectionValues$N_11,
                                 infectionValues$N_01,
                                 infectionValues$N_cntrl_group))
  
  colnames(weights) <- c("nseed", "sum1", "sum2","sum3","sum4")
  # take the sums needed
  sums1 <- aggregate(.~nseed,weights, function(x) sum = sum(x))
  # repopulate with a sums column the same length as the df
  sums <- merge(sums1, infectionValues, by = "nseed")
  
  #infectionValues$sum1 <- tapply(infectionValues$N_trt_group,INDEX = infectionValues$nseed, sum)
  
  
  # two stage inverse probability weights
  # use sums from above
  sums$w11 <- (sums$sum1/N)^(-1) *(sums$sum2/sums$sum1)^(-1)
  #assign("test", sums$w11,envir = .GlobalEnv)
  
  sums$w01 <- (sums$sum1/N)^(-1) * (sums$sum3/sums$sum1)^(-1)
  sums$w0 <- (sums$sum4/N)^(-1)
  sums$w1 <- (sums$sum1/N)^(-1)
  
  #make sure the weights are in the correct order for data set
  sums2<-as.data.frame(cbind(sums$nseed,sums$w11,sums$w01,sums$w0,sums$w1))
  names(sums2) <- c("nseed", "w11","w01","w0","w1")
  #keep one record per nseed
  sums3 <- sums2 %>% 
    group_by(nseed) %>% 
    filter(row_number()==1) %>%
    as.data.frame()
  #sums2a <- do.call(rbind, lapply(split(sums2, sums2$nseed), head, 1))
  
  
  #assign("wcheck", w1, envir = .GlobalEnv)
  Y_prev <- infectionValues$HIV_11 * infectionValues$N_11 * infectionValues$N_trt_group +
    infectionValues$HIV_01 * (1-infectionValues$N_11) *infectionValues$N_trt_group +
    infectionValues$HIV_00 * (1-infectionValues$N_11) * (1-infectionValues$N_trt_group)
  
  Y_inc <- infectionValues$Inc_11 * infectionValues$N_11 * infectionValues$N_trt_groupi +
    infectionValues$Inc_01 * (1-infectionValues$N_11) *infectionValues$N_trt_groupi +
    infectionValues$Inc_00 * (1-infectionValues$N_11) * (1-infectionValues$N_trt_groupi)

  infectionValues2 <- merge(sums3,infectionValues, by = "nseed")
  infectionValues3 <- merge(sums3i,infectionValues, by = "nseed")
 
  #risks
  # this calculates them separately then finds the sum
  risk11_prev <- infectionValues2$w11*infectionValues2$w_h*infectionValues2$HIV_11 #*(infectionValues$N_11*infectionValues$N_trt_group)
  risk11_prev <- tapply(risk11_prev, infectionValues2$nseed, sum)
  #Some of these estimates are out of bounds >1?
  risk01_prev <- infectionValues2$w01*infectionValues2$w_h*infectionValues2$HIV_01#*(infectionValues$N_01)
  risk01_prev <- tapply(risk01_prev, infectionValues2$nseed, sum)
  #Some of these estimates are out of bounds >1?
  risk00_prev <- infectionValues2$w0*infectionValues2$w_h*infectionValues2$HIV_00 #*(1-infectionValues$N_11)*(1 - infectionValues$N_trt_group)
  risk00_prev <- tapply(risk00_prev, infectionValues2$nseed, sum)
  # 
  risk1_prev <- infectionValues2$w1 * infectionValues2$w_h * (infectionValues2$HIV_11 + infectionValues2$HIV_01) #* infectionValues$N_trt_group
  risk1_prev <- tapply(risk1_prev, infectionValues2$nseed, sum)
  risk0_prev <- risk00_prev
  #incidence risk
  risk11_inc <- infectionValues3$w11*infectionValues3$w_h*infectionValues3$Inc_11 #*infectionValues$N_11i*infectionValues$N_trt_groupi
  risk11_inc <- tapply(risk11_inc, infectionValues3$nseed, sum)
  
  risk01_inc <- infectionValues3$w01*infectionValues3$w_h*infectionValues3$Inc_01 #*(infectionValues$N_01i)
  risk01_inc <- tapply(risk01_inc, infectionValues3$nseed, sum)
  
  risk00_inc <- infectionValues3$w0*infectionValues3$w_h*infectionValues3$Inc_00 #*(1-infectionValues$N_11)*(1 - infectionValues$N_trt_group)
  risk00_inc <- tapply(risk00_inc, infectionValues3$nseed, sum)
  # 
  risk1_inc <- infectionValues3$w1 * infectionValues3$w_h * (infectionValues3$Inc_11 + infectionValues3$Inc_01) #* infectionValues$N_trt_groupi
  risk1_inc <- tapply(risk1_inc, infectionValues3$nseed, sum)
  risk0_inc <- risk00_inc
  

  risk <- as.data.frame(cbind(risk11_prev,risk01_prev,risk00_prev,risk1_prev,risk0_prev, risk11_inc,risk01_inc,risk00_inc,risk1_inc,risk0_inc))
  title <- paste("risk", coverage, sep = "")
  
  assign(title, risk, envir = .GlobalEnv)
}

#Risks must be between 0 and 1 - with weights are some out of bounds?
#riskDF_10[ which(riskDF_10$risk01_prev>1 ), ]

#infectionValues[which(infectionValues$nseed=="9317"),]


riskDF_10<-estimators(treatment10,10)
riskDF_20<-estimators(treatment20,20)
riskDF_30<-estimators(treatment30,30)
riskDF_40<-estimators(treatment40,40)
riskDF_50<-estimators(treatment50,50)
riskDF_60<-estimators(treatment60,60)
riskDF_70<-estimators(treatment70,70)
riskDF_80<-estimators(treatment80,80)
riskDF_90<-estimators(treatment90,90)

risk10ave <- lapply(riskDF_10,2,FUN=mean)
risk20ave <- lapply(riskDF_20,2,FUN=mean)
risk30ave <- lapply(riskDF_30,2,FUN=mean)
risk40ave <- lapply(riskDF_40,2,FUN=mean)
risk50ave <- lapply(riskDF_50,2,FUN=mean)
risk60ave <- lapply(riskDF_60,2,FUN=mean)
risk70ave <- lapply(riskDF_70,2,FUN=mean)
risk80ave <- lapply(riskDF_80,2,FUN=mean)
risk90ave <- lapply(riskDF_90,2,FUN=mean)

write.table(risk10ave,"risk10.txt",row.names = F)
write.table(risk20ave,"risk20.txt",row.names = F)
write.table(risk30ave,"risk30.txt",row.names = F)
write.table(risk40ave,"risk40.txt",row.names = F)
write.table(risk50ave,"risk50.txt",row.names = F)
write.table(risk60ave,"risk60.txt",row.names = F)
write.table(risk70ave,"risk70.txt",row.names = F)
write.table(risk80ave,"risk80.txt",row.names = F)
write.table(risk90ave,"risk90.txt",row.names = F)


