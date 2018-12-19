N <- 11245

estimators <- function(treatment,coverage){
  treatment$treatmentgp <- (treatment$Nprep + treatment$NprepElig)/(treatment$Nprep + treatment$NprepElig)
  treatment[is.na(treatment)] <-0
  time0 <- subset(treatment,t == 0)
  time60 <- subset(treatment, t == max(treatment$t))

  
  #start_01 <- subset(time0, NprepElig != 0) # number of people not on PrEP in treatment gp
  #
  #end_01 <- subset(time60, time0$NprepElig != 0)
  N_01 <- time0$NprepElig
  HIV_01 <- time0$NprepElig - time60$NprepElig + time0$Nhiv*time0$treatmentgp # no. of people who
          # lose prep elig = no. prep elig infected (plus those who started HIV+)
  assign("test", HIV_01, envir = .GlobalEnv)
  Inc_01 <- time0$NprepElig - time60$NprepElig # no. of people who
  # lose prep elig = no. prep elig infected
  N_11 <- time0$Nprep # number of people on PrEP
  HIV_11 <- time60$Nhiv * time60$treatmentgp - HIV_01
  Inc_11 <- time60$Nnewinf * time60$treatmentgp - Inc_01
  
  N_trt_group <- N_11 + N_01
  
  #N_00 <- time0$N[treatment$NprepElig == 0 & treatment$Nprep == 0 & treatment$t == 0] # population of control gp
  HIV_00 <- time60$Nhiv - HIV_11 - HIV_01 # HIV population of control gp
  Inc_00 <- time60$Nnewinf - Inc_11 - Inc_01
  
  #individual weights
  K <- with(time0,ave(compID,nseed,FUN = length))
  
  w_h <-  1/(time0$totalN*K)
  w_i <- 1/N
  
  #two stage inverse probability weights
  w11 <-(sum(N_trt_group)/N)^(-1)*(sum(N_11*N_trt_group)/sum(N_trt_group))^(-1)
  w01 <-(sum(N_trt_group)/N)^(-1)*(sum((1 - N_11)*N_trt_group)/sum(N_trt_group))^(-1)
  w0 <-(sum((1 - N_trt_group))/N)^(-1)
  w1 <- (sum(N_trt_group)/N)^(-1)
  
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

