temp2 <- read.table(temp[1], header = T)
datalist[[1]] <- temp2
temp2 <- read.table(temp[2], header = T)
datalist[[2]] <- temp2
temp2 <- read.table(temp[3], header = T)
datalist[[3]] <- temp2
temp2 <- read.table(temp[4], header = T)
datalist[[4]] <- temp2
temp2 <- read.table(temp[5], header = T)
datalist[[5]] <- temp2
temp2 <- read.table(temp[6], header = T)
datalist[[6]] <- temp2
temp2 <- read.table(temp[7], header = T)
datalist[[7]] <- temp2
temp2 <- read.table(temp[8], header = T)
datalist[[8]] <- temp2
temp2 <- read.table(temp[9], header = T)
datalist[[9]] <- temp2
temp2 <- read.table(temp[10], header = T)
datalist[[10]] <- temp2


for (i in 1:9){
  y <- paste("t2_", i, "0", sep = "")
  assign(x =y, as.data.frame(datalist[i]), envir = .GlobalEnv)}
