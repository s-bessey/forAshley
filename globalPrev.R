t0 <- subset(treatment50, treatment50$t ==0)
t0$prepgp <- ifelse(t0$Nprep + t0$NprepElig > 0, 1,0)
t60 <- subset(treatment50, treatment50$t > 0)

test50 <- treatment50
test50$prepgp <-NA
test50$prepgp <- ifelse(test50$t == 0 & test50$Nprep+test50$NprepElig >0, 1, test50$prepgp)
test50$prepgp <- ifelse(test50$t == 0 & test50$Nprep+test50$NprepElig== 0, 0, test50$prepgp)
test50 <- test50[order(test50$compID),]
test50 <- test50[order(test50$nseed, test50$compID),]
test50$prepgp <- na.locf(test501$prepgp)

globalPrev <- function(treatment){
  #need to get sum of each used variable by seed
  # prep <- subset(treatment, treatment$prepgp > 0)
  #prep0 <- subset(treatment, treatment$t == 0)
  # prepN <- sum(prep0$totalN)
  # prep0 <- aggregate(.~nseed, prep0, function(x) sum = sum(x))
  # prepHIV0 <- mean(prep0$Nhiv)
  prepOutput <- aggregate(.~nseed + t + prepgp, treatment, FUN = sum)
  prepOutput <- aggregate(.~t + prepgp, prepOutput, function(x) mean=mean(x))
  #prepOutput <- mean(treatment$totalN[treatment$t ==0])
  assign("test", prepOutput, envir = .GlobalEnv)
  prep <- subset(treatment, treatment$Nprep + treatment$NprepElig == 0)
  prepOutput2 <- aggregate(.~nseed + t, prep, FUN = sum)
  prepOutput2 <- aggregate(.~t, prepOutput2, FUN = mean)
  assign("test2", prepOutput2, envir = .GlobalEnv)
  
}

  
  prep60 <- subset(prep, prep$t == 60)
  prep60 <- aggregate(prep60$nseed, prep60, function(x) sum = sum(x))
  prep60$t <- rep(60, nrow(prep60))
  prep60 <- aggregate(.~t, prep60, function(x) mean = mean(x))
  prepHIV60 <- mean(prep60$Nhiv)

  noprep <- subset(treatment, treatment$prepgp == 0)
  noprep0 <- subset(prep, prep$t == 0)
  noprepN <- sum(noprep0$totalN)
  noprep0 <- aggregate(.~nseed, noprep0, function(x) sum = sum(x))
  noprep0 <- aggregate(.~t, noprep0, function(x) mean = mean(x))
  noprep60 <- subset(noprep, noprep$t == 60)
  noprep60 <- aggregate(.~nseed, noprep60, function(x) sum = sum(x))
  noprep60$t <- rep(60, nrow(noprep60))
  noprep60 <- aggregate(.~t, noprep60, function(x) mean = mean(x))


  t <- prep$t
  prevPrep <- c(prep0$Nhiv/prepN, prep60$Nhiv/prepN)
  incPrep <- c((prep60$Nhiv-prep0$Nhiv)/prepN)
  prevNoPrep <- c(noprep0$Nhiv/noprepN, noprep60$Nhiv/noprepN)
  incNoPrep <- c((noprep60$Nhiv-noprep0$Nhiv)/noprepN)

  assign("prevalence", as.data.frame(cbind(prevPrep, incPrep, prevNoPrep,incNoPrep)), envir = .GlobalEnv)
  assign("nhiv", c(prepHIV0, prepHIV60, prepN, noprep0$totalN), envir = .GlobalEnv)
  return(prep60)
}

globalPrev(allTime2)
fig <- ggplot(prevalence, aes(x = t, y = prevPrep)) +
  geom_line(size = 1) + geom_line(aes(x = t, y = prevNoPrep), size = 1) +
  xlab("t")+ ylab("Prevalence (%)") + 
  theme_classic()

noprepID <- subset(treatment50, treatment50$NprepElig+treatment50$Nprep == 0)
noprepID <- unique(noprepID$compID)

noprep <- subset(treatment50, treatment50$NprepElig+treatment50$Nprep == 0)
prep <- subset(treatment50, treatment50$NprepElig+treatment50$Nprep != 0)

aggtNoPrep <- aggregate(.~nseed+t, noprep, function(x) sum(x))
timesNoPrep <- aggregate(.~t, aggtNoPrep, FUN = mean)
prevFigNoPrep <- ggplot(timesNoPrep) + geom_line(aes(x = timesNoPrep$t, y = timesNoPrep$Nhiv/timesNoPrep$totalN))


incFig <- ggplot(times) + geom_line(aes(x=times$t, y = times$Nnewinf/times$totalN))

