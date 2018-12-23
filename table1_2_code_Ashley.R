


treatmentCoverages <- as.data.frame(c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%"))


table13function <- function(treatment, coverages){

  treatment$treatmentgp <- (treatment$Nprep + treatment$NprepElig)/(treatment$Nprep + treatment$NprepElig)
  treatment[is.na(treatment)] <- 0
  

  time0 <- subset(treatment, t ==0)
  time60 <- subset(treatment, t != 0)

  
  time60$HIV_01 <- time0$NprepElig - time60$NprepElig + (time0$Nhiv * time0$treatmentgp)
  time60$HIV_11 <- (time60$Nhiv * time0$treatmentgp) - time60$HIV_01
  time60$inc_01 <- time0$NprepElig - time60$NprepElig
  time60$inc_11 <- (time60$Nnewinf * time0$treatmentgp) - time60$inc_01
  time60$N_11 <- time0$Nprep
  time60$N_01 <- time0$NprepElig
  
  total_Ns <- aggregate(.~nseed, time60, sum)
  N_11 <- mean(total_Ns$N_11)
  N_01 <- mean(total_Ns$N_01)
  HIV_01 <- mean(total_Ns$HIV_01)
  HIV_11 <- mean(total_Ns$HIV_11)
  Inc_01 <- mean(total_Ns$inc_01)
  Inc_11 <- mean(total_Ns$inc_11)
  
  t1 <- as.data.frame(cbind(coverages, N_11, HIV_11, HIV_11/N_11, N_01, HIV_01, HIV_01/N_01))
  colnames(t1) <- c("Component PrEP Coverage Level", "Total Persons", "HIV+", "Point Prevalence, %", "Total Persons", "HIV+", "Point Prevalence, %")
  t3 <- cbind(coverages, N_11, Inc_11, Inc_11/N_11, N_01, Inc_01, Inc_01/N_01)
  name <- paste("t1",coverages,sep = "_")
  assign(name, as.data.frame(t1), envir = .GlobalEnv)
  name <- paste("t2",coverages,sep = "_")
  assign(name, as.data.frame(t2), envir = .GlobalEnv)
}


table13function(treatment10, 10)
table13function(treatment20, 20)
table13function(treatment30, 30)
table13function(treatment40, 40)
table13function(treatment50, 50)
table13function(treatment60, 60)
table13function(treatment70, 70)
table13function(treatment80, 80)
table13function(treatment90, 90)


table1 <- as.data.frame(rbind(t1_10, t1_20, t1_30, t1_40, t1_50, t1_60, t1_70, t1_80, t1_90))
# table1 <- as.data.frame(cbind(treatmentCoverages, table1))
tab_df(table1, col.header = colnames(table1),
       show.rownames = F, alternate.rows = T,file="Table1.doc")

table3 <- as.data.frame(rbind(t2_10, t2_20, t2_30, t2_40, t2_50, t2_60, t2_70, t2_80, t2_90))

tab_df(table3, col.header = colnames(table1),
       show.rownames = F, alternate.rows = T,file="Table3.doc")
