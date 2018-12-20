N <- 11245
setwd("Ashley")
temp=list.files(pattern = "*10*")
#treatment10 <- data.frame(rseed = integer(), pseed = integer(),nseed = integer(),
 #                         t = integer(), compID = integer(), totalN = integer(),
  #                        Nhiv = integer(), NprepElig = integer(), Nprep = integer(),
   #                       Nnewinf = integer())
datalist <- list()
temp=list.files(pattern = "*10_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(temp[i],header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment10 <- do.call(rbind, datalist)

temp=list.files(pattern = "*20_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment20 <-docall(rbind,datalist)

temp=list.files(pattern = "*30_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment30 <-docall(rbind,datalist)

temp=list.files(pattern = "*40_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment40 <-docall(rbind,datalist)

temp=list.files(pattern = "*50_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment50 <-docall(rbind,datalist)

temp=list.files(pattern = "*60_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment60 <-docall(rbind,datalist)

temp=list.files(pattern = "*70_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment70 <-docall(rbind,datalist)

temp=list.files(pattern = "*80_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment80 <-docall(rbind,datalist)

temp=list.files(pattern = "*90_*")
for (i in 1:length(temp)){
  assign(paste("temp2",sep = ""), read.table(location,header = TRUE))
  #assign(temp[i],temp2)
  datalist[[i]] <- temp2
}
treatment90 <-docall(rbind,datalist)








estimators <- function(treatment,coverage){
  treatment$treatmentgp <- (treatment$Nprep + treatment$NprepElig)/(treatment$Nprep + treatment$NprepElig)
  treatment[is.na(treatment)] <-0
  time0 <- subset(treatment,t == 0)
  time60 <- subset(treatment, t == max(treatment$t))
  infectionValues <- as.data.frame(time60$nseed)
  infectionValues$compID <- time60$compID
  
  #start_01 <- subset(time0, NprepElig != 0) # number of people not on PrEP in treatment gp
  #
  #end_01 <- subset(time60, time0$NprepElig != 0)
  infectionValues$N_01 <- time0$NprepElig
  infectionValues$HIV_01 <- time0$NprepElig - time60$NprepElig + time0$Nhiv*time0$treatmentgp # no. of people who
          # lose prep elig = no. prep elig infected (plus those who started HIV+)
  assign("test", infectionValues$HIV_01, envir = .GlobalEnv)
  infectionValues$Inc_01 <- time0$NprepElig - time60$NprepElig # no. of people who
  # lose prep elig = no. prep elig infected
  infectionValues$N_11 <- time0$Nprep # number of people on PrEP
  infectionValues$HIV_11 <- time60$Nhiv * time60$treatmentgp - infectionValues$HIV_01
  infectionValues$Inc_11 <- time60$Nnewinf * time60$treatmentgp - infectionValues$Inc_01
  
  infectionValues$N_trt_group <- infectionValues$N_11 + infectionValues$N_01
  
  #N_00 <- time0$N[treatment$NprepElig == 0 & treatment$Nprep == 0 & treatment$t == 0] # population of control gp
  infectionValues$HIV_00 <- time60$Nhiv - infectionValues$HIV_11 - infectionValues$HIV_01 # HIV population of control gp
  infectionValues$Inc_00 <- time60$Nnewinf - infectionValues$Inc_11 - infectionValues$Inc_01
  
  #individual weights
  K <- with(time0,ave(compID,nseed,FUN = length))
  
  infectionValues$w_h <-  1/(time0$totalN*K)
  infectionValues$w_i <- 1/N
  
  #two stage inverse probability weights
  w11 <-(sum(N_trt_group)/N)^(-1)*(sum(N_11*N_trt_group)/sum(N_trt_group))^(-1)
  #sum(number of people in treatment groups)/N(^-1) *sum(nprep*N_trt_group)
  w01 <-(sum(N_trt_group)/N)^(-1)*(sum((1 - N_11)*N_trt_group)/sum(N_trt_group))^(-1)
  w0 <-(sum((1 - N_trt_group))/N)^(-1)
  w1 <- (sum(N_trt_group)/N)^(-1)
  
  #functions for risks
  
  w11fun <- function(weightInput){
    weight11<-(sum(weightInput$N_trt_gp)/N)^(-1)*(sum(weightInput$N_11*weightInput$N_trt_group)/sum(weightInput$trtgp))^(-1)
    assign("w11",weight11)
  }

  #risks
  risk11_prev <-sum(w11*w_h*HIV_11*N_11*N_trt_group) 
  risk01_prev <-sum(w01*w_h*HIV_01*(1-N_11)*N_trt_group)
  risk00_prev <-sum(w0*w_h*HIV_00*(1-N_11)*(1 - N_trt_group))
  #risk00_prev <- summarize(group_by(nseed,compID),sum)
  risk1_prev <-sum(w1 * w_h * (HIV_11 + HIV_01) * N_trt_group)
  risk0_prev <-sum(w0*w_h*HIV_00*(1 - N_trt_group))
  
  risk11_inc <-sum(w11 * w_h * Inc_11 * N_11 * N_trt_group)
  risk01_inc <-sum(w01 * w_h * Inc_01 * (1-N_11) * N_trt_group)
  risk00_inc <-sum(w0 * w_h * Inc_00 * (1-N_11) * (1 - N_trt_group))
  
  risk1_inc <-sum(w1 * w_h * (Inc_11+Inc_01) * N_trt_group)
  risk0_inc <-sum(w0 * w_h * Inc_00 * (1 - N_trt_group))
  
  risk <- cbind(risk11_prev,risk01_prev,risk00_prev,risk1_prev,risk0_prev, risk11_inc,risk01_inc,risk00_inc,risk1_inc,risk0_inc)
  title <- paste("risk", coverage, sep = "")
  assign(title,risk, envir = .GlobalEnv)
}

