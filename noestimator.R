treatment10 <- as.data.frame(rbind(read.table("Ashley/10_10Report.txt",header = T), read.table("Ashley/10_9Report.txt",header = T),
                     read.table("Ashley/10_8Report.txt",header = T), read.table("Ashley/10_7Report.txt", header = T)))
tester <- function(treatment,coverage){
  treatment$treatmentgp <- (treatment$Nprep + treatment$NprepElig)/(treatment$Nprep + treatment$NprepElig)
  treatment[is.na(treatment)] <-0
  time0 <- subset(treatment,t == 0)
  time60 <- subset(treatment, t == max(treatment$t))
  
  
  #start_01 <- subset(time0, NprepElig != 0) # number of people not on PrEP in treatment gp
  #
  #end_01 <- subset(time60, time0$NprepElig != 0)
  N_01 <- time0$NprepElig + (time0$Nhiv * time0$treatmentgp)
  HIV_01 <- (time0$NprepElig - time60$NprepElig) + (time0$Nhiv * time0$treatmentgp) # no. of people who
  # lose prep elig = no. prep elig infected (plus those who started HIV+)
  assign("test", HIV_01, envir = .GlobalEnv)
  Inc_01 <- time0$NprepElig - time60$NprepElig # no. of people who
  # lose prep elig = no. prep elig infected
  N_11 <- time0$Nprep # number of people on PrEP
  HIV_11 <- (time60$Nhiv * time0$treatmentgp) - HIV_01
  Inc_11 <- time60$Nnewinf * time0$treatmentgp - Inc_01
  
  N_trt_group <- N_11 + N_01
  
  #N_00 <- time0$N[treatment$NprepElig == 0 & treatment$Nprep == 0 & treatment$t == 0] # population of control gp
  HIV_00 <- time60$Nhiv - HIV_11 - HIV_01 # HIV population of control gp
  Inc_00 <- time60$Nnewinf - Inc_11 - Inc_01
  N_00 <- time0$totalN * (1-time0$treatmentgp)
  
  risk_11 <- HIV_11/N_11
  risk_11[is.na(risk_11)] <-0
  risk_01 <- HIV_01/N_01
  risk_01[is.na(risk_01)] <-0
  risk_00 <- HIV_00/N_00
  risk_1 <- (HIV_11+HIV_01)/(N_01+N_11)
  assign("risk10", as.data.frame(cbind(risk_11,risk_01,risk_1,risk_00)),envir = .GlobalEnv)
  
  }
