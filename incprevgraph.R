treatment$treatmentgp <- (treatment$Nprep + treatment$NprepElig)/(treatment$Nprep + treatment$NprepElig)
treatment[is.na(treatment)] <-0

# issues: on/off prep? no, because ummm only need it at t = 0
N_01_prev <- totalN - Nhiv
N_01_inc <- NprepElig[t == 0]
HIV_01
cumInf_01
N_00 <- N
HIV_00
cumInf_00


incidenceCalc <- function(treatment){
  treatment$treatmentgp <- (treatment$Nprep + treatment$NprepElig)/(treatment$Nprep + treatment$NprepElig)
  treatment[is.na(treatment)] <-0
  treatment$prevalence <- treatment$Nhiv / treatment$totalN
  
  A0 <- subset(treatment,treatment$treatmentgp == 0)
  A1 <- subset(treatment, treatment$treatmentgp != 0)

  x <-treatment %>% group_by(nseed) %>% mutate(cumInf = cumsum(Incid))
  forCum <- data.frame(x)
}